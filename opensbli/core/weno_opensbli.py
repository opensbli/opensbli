from sympy import IndexedBase, Symbol, pprint, Rational, solve, interpolating_poly, integrate, Eq, sqrt, zeros, Abs, Float, Matrix, flatten, Max, diag, Function
from opensbli.core.opensblifunctions import WenoDerivative
from opensbli.core.opensbliobjects import EinsteinTerm, DataSetBase, ConstantObject
from opensbli.core.opensbliequations import SimulationEquations
from opensbli.core.kernel import Kernel
from opensbli.core.grid import GridVariable
from opensbli.utilities.helperfunctions import increment_dataset as incr_dset
from opensbli.physical_models.ns_physics import NSphysics
from .scheme import Scheme


class WenoHalos(object):
    def __init__(self, order, reconstruction=None):
        """ Object for WENO halos.
        arg: int: order: Order of the WENO scheme.
        arg: bool: reconstruction: True if halos for a reconstruction. """
        k = int(0.5*(order+1))
        if not reconstruction:
            self.halos = [-k, k+1]
        else:
            self.halos = [-1, 1]
        return

    def get_halos(self, side):
        return self.halos[side]


class ConfigureWeno(object):
    def __init__(self, k, side):
        """ Object containing the parameters needed by the WENO reconstruction for a given side.
        arg: int: k: WENO coefficient k, equal to scheme order = 2k - 1.
        arg: int: side: Side of the WENO reconstruction, either -1 (left) or +1 (right)."""
        self.side = side
        if self.side == -1:
            self.name, self.short_name, self.shift = 'left', 'L', 1
        elif self.side == 1:
            self.name, self.short_name, self.shift = 'right', 'R', 0
        self.k, self.side = k, side
        self.func_points = self.generate_left_right_points()
        # k passed explicitly as 2 sets of ENO coefficients are needed for the smoothness indicators
        self.c_rj = self.generate_eno_coefficients(k)
        self.c2_rj = self.generate_eno_coefficients(2*k-1)
        self.opt_weights = self.generate_optimal_weights()
        print self.opt_weights
        return

    @property
    def generate_symbolic_function_points(self):
        """ Creates a set of symbolic function points f[0], f[1], ... to be used in the
        WENO reconstruction, and a dictionary linking each function point to its index.
        :returns list symbolic_functions: The symbolic function points.
        :returns dict symbolic_function_dictionary: A dictionary of index:symbolic function pairs."""
        f = IndexedBase('f')
        symbolic_functions = []
        symbolic_function_dictionary = {}
        for p in set(self.func_points):
            symbolic_functions.append(f[p])
            symbolic_function_dictionary[p] = symbolic_functions[-1]
        return symbolic_functions, symbolic_function_dictionary

    def generate_eno_coefficients(self, k):
        """ Generates the c_rj ENO reconstruction coefficients given in Table 2.1
        of 'Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes
        for Hyperbolic Conservation Laws' by Shu(1997). Computation is of equation (2.21).

        :arg int: k: WENO coefficient k, equal to scheme order = 2k - 1.
        :returns: dict: c_rj: Dictionary in the form (r,j) : ENO coefficient c_rj."""
        if self.side == 1:
            d = 1
        elif self.side == -1:
            d = 0
        r_values, j_values = range(k), range(k)
        c_rj = {}
        for r in r_values:
            for j in j_values:
                c_rj_sum = 0
                for m in range(j+1, k+1):
                    top_sum = 0
                    bottom_product = 1
                    for l in [x for x in range(k+1) if x != m]:
                        top_product_total = 1
                        for q in [x for x in range(k+1) if (x != m and x != l)]:
                            top_product = r - q + d
                            top_product_total *= top_product
                        bottom_product *= m - l
                        top_sum += top_product_total
                    c_rj_sum += Rational(top_sum, bottom_product)
                c_rj[(r, j)] = c_rj_sum
        return c_rj

    def generate_optimal_weights(self):
        """ Creates the optimal weights d_r needed by the WENO scheme, as specified in Shu(1997).
        :returns dict: coeff_dict: A dictionary of the optimal weights. """
        # Store reconstruction coefficients to form the sum
        k, c_rj, c2_rj = self.k, self.c_rj, self.c2_rj
        opt_weights = [Symbol('d_%d' % r) for r in range(k)]
        r_values = [[k-1-j for j in range(m+1)][::-1] for m in range(k)]
        equations = []
        for i, indices in enumerate(r_values):
            variables = [opt_weights[r] for r in indices]
            coeffs = [c_rj[(r, j)] for j, r in enumerate(indices)]
            rhs = c2_rj[(k-1, i)]
            v = [variables[i]*coeffs[i] for i in range(len(indices))]
            equations += [sum(v) - rhs]
        # Solve the linear system to get the coefficients
        solution = solve(equations, opt_weights)
        coeff_dict = {}
        for i in range(k):
            coeff_dict[(0, i)] = solution[opt_weights[i]]
        return coeff_dict

    def generate_left_right_points(self):
        """ Evaluate the required function evaluation points for left and
        right WENO reconstructions.
        returns: list: fn_points: A list of integers to index the function locations."""
        k = self.k
        if self.side == -1:
            c = 1
        elif self.side == 1:
            c = 0
        fn_points = []
        for r in range(k):
            fn_points += [-r+j+c for j in range(k)]
        return fn_points


class JS_smoothness(object):
    def __init__(self, k):
        """ An object to hold the Jiang-Shu smoothness coefficients. Listed as Eq(3.1) in
        'WENO coefficient k, equal to scheme order = 2k - 1.' Shu(1996)."""
        self.k = k
        self.smooth_coeffs = self.generate_smoothness_coefficients()
        return

    def generate_smoothness_coefficients(self):
        """Extracts the JS smoothness coefficients.
        returns: dict: Dictionary of coefficients used in the generation of the smoothness indicators
        calculated in the 'generate_function_smoothness_indicators' routine."""
        k = self.k
        smooth_coeffs = {}
        x, dx = Symbol('x'), Symbol('dx')
        for r in range(k):
            # Generate x, y values to create interpolating polynomial
            dx_values = [dx*i for i in range(-r, k+1-r)]
            nodes = [Symbol('fn[i%+d]' % i) for i in range(-r, k-r)]
            funcs = [0]
            for i in range(k):
                funcs.append(funcs[-1] + dx*nodes[i])
            # Perform Lagrange interpolation
            lagrange_interp = interpolating_poly(k+1, x, dx_values, funcs).diff(x)
            # Loop over the derivatives of the interpolating polynomial
            total = 0
            for l in range(1, k):
                q = (lagrange_interp.diff(x, l))**2
                # Perform the integration and multiply by h^(2*l-1) over cell
                q = integrate(q.as_poly(x), x) * dx**(2*l-1)
                total += (q.subs(x, dx) - q.subs(x, 0))
            done = []
            # Save the coefficients of the smoothness indicator
            for m in range(0, 2*k-2):
                for n in range(0, 2*k-2):
                    func_product = Symbol('fn[i%+d]' % (-r+m))*Symbol('fn[i%+d]' % (-r+n))
                    if func_product not in done:
                        c = total.coeff(func_product)
                    if c != 0:
                        smooth_coeffs[(r, m, n)] = c
                    done.append(func_product)
        return smooth_coeffs

    def generate_function_smoothness_indicators(self, fn):
        """ Generates the standard Jiang-Shu smoothness indicators as polynomials in the
        symbolic placeholder functions 'f'.
        arg: object: fn: The reconstruction variable object to contain the smoothness indicators."""
        if isinstance(fn, LeftWenoReconstructionVariable):
            shift = 1
        elif isinstance(fn, RightWenoReconstructionVariable):
            shift = 0
        k = self.k
        # Compute the smoothness indicator and alpha
        for r in range(k):
            fn.smoothness_symbols += [Symbol('beta_%d' % r)]
            local_smoothness = 0
            for m in range(0, 2*k-2):  # WARNING: change to (0, 2*k -1) if this is now broken.
                for n in range(0, 2*k-2):
                    beta = self.smooth_coeffs.get((r, m, n))
                    if beta is not None:
                        shift_indices = [-r+m+shift, -r+n+shift]
                        func_product = fn.function_stencil_dictionary[shift_indices[0]]*fn.function_stencil_dictionary[shift_indices[1]]
                        local_smoothness += beta*func_product
            pprint(local_smoothness)
            # print "old: ", local_smoothness.count_ops()
            # local_smoothness = horner(local_smoothness)
            # print "new: ", local_smoothness.count_ops()
            # pprint(local_smoothness)
            # print "\n\n"
            fn.smoothness_indicators += [local_smoothness]
        return


class WenoReconstructionVariable(object):
    def __init__(self, name):
        """ Reconstruction variable object to hold the quantities required for WENO.
        arg: str: name: Name of the reconstruction, either left or right."""
        self.name = name
        self.smoothness_indicators = []
        self.smoothness_symbols = []
        self.alpha_evaluated = []
        self.alpha_symbols = []
        self.omega_evaluated = []
        self.omega_symbols = []
        self.function_stencil_dictionary = {}
        self.reconstructed_expression = None
        self.reconstructed_symbol = GridVariable('%s_%s' % ('reconstruct', name))
        return

    def update_quantities(self, original):
        """ Updates the quantities required by WENO in the reconstruction variable.
        arg: object: original: Reconstruction object variable, either left or right reconstruction."""
        print "Reconstruction variable: ", self.name
        self.smoothness_symbols += [GridVariable('%s%s' % (s, self.name)) for s in original.smoothness_symbols]
        self.alpha_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.alpha_symbols]
        self.omega_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.omega_symbols]

        subs_dict = dict(zip(original.smoothness_symbols+original.alpha_symbols+original.omega_symbols, self.smoothness_symbols+self.alpha_symbols+self.omega_symbols))
        for key, value in original.function_stencil_dictionary.iteritems():
            subs_dict[value] = self.function_stencil_dictionary[key]

        self.smoothness_indicators = [s.subs(subs_dict) for s in original.smoothness_indicators]
        self.alpha_evaluated = [s.subs(subs_dict) for s in original.alpha_evaluated]
        self.omega_evaluated = [s.subs(subs_dict) for s in original.omega_evaluated]
        self.reconstructed_expression = original.reconstructed_expression.subs(subs_dict)
        return

    def add_evaluations_to_kernel(self, kernel):
        all_symbols = self.smoothness_symbols + self.alpha_symbols + self.omega_symbols
        all_evaluations = self.smoothness_indicators + self.alpha_evaluated + self.omega_evaluated
        for no, value in enumerate(all_symbols):
            kernel.add_equation(Eq(value, all_evaluations[no]))
        kernel.add_equation(Eq(self.reconstructed_symbol, self.reconstructed_expression))
        return


class LeftWenoReconstructionVariable(WenoReconstructionVariable):
    def __init__(self, name):
        """ Reconstruction object for the left reconstruction.
        arg: str: name: 'left' """
        WenoReconstructionVariable.__init__(self, name)
        return


class RightWenoReconstructionVariable(WenoReconstructionVariable):
    def __init__(self, name):
        """ Reconstruction object for the right reconstruction.
        arg: str: name: 'right' """
        WenoReconstructionVariable.__init__(self, name)
        return


class WenoZ(object):
    def __init__(self, k):
        self.k = k
        self.eps = 1e-16
        return

    def generate_alphas(self, RV, WenoConfig):
        """ Create the alpha terms for the non-linear WENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        tau_5 = Abs(RV.smoothness_symbols[0] - RV.smoothness_symbols[-1])
        p = 2
        for r in range(self.k):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append(WenoConfig.opt_weights.get((0, r))*(Float(1) + (tau_5/(self.eps + RV.smoothness_symbols[r]))**p))
        return

    def generate_omegas(self, RV, WenoConfig):
        """ Create the omega terms for the non-linear WENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        for r in range(self.k):
            RV.omega_symbols += [Symbol('omega_%d' % r)]
            RV.omega_evaluated += [RV.alpha_symbols[r]/sum(RV.alpha_symbols)]
        return

    def generate_reconstruction(self, RV, WenoConfig):
        """ Create the final WENO stencil by summing the stencil points, ENO coefficients and WENO weights..
        arg: object RV: The reconstruction variable object.
        arg: object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        stencil = 0
        k, c_rj = self.k, WenoConfig.c_rj
        for r in range(k):
            eno_expansion = [c_rj.get((r, j))*RV.function_stencil_dictionary[WenoConfig.func_points[r*k+j]] for j in range(k)]
            stencil += RV.omega_symbols[r]*sum(eno_expansion)
        RV.reconstructed_expression = stencil
        return


class WenoJS(object):
    def __init__(self, k):
        self.k = k
        self.eps = 1e-6
        return

    def generate_alphas(self, RV, WenoConfig):
        """ Create the alpha terms for the non-linear WENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        for r in range(self.k):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append(WenoConfig.opt_weights.get((0, r))/(self.eps+RV.smoothness_symbols[r])**2)
        return

    def generate_omegas(self, RV, WenoConfig):
        """ Create the omega terms for the non-linear WENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        for r in range(self.k):
            RV.omega_symbols += [Symbol('omega_%d' % r)]
            RV.omega_evaluated += [RV.alpha_symbols[r]/sum(RV.alpha_symbols)]
        return

    def generate_reconstruction(self, RV, WenoConfig):
        """ Create the final WENO stencil by summing the stencil points, ENO coefficients and WENO weights..
        arg: object RV: The reconstruction variable object.
        arg: object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        stencil = 0
        k, c_rj = self.k, WenoConfig.c_rj
        for r in range(k):
            eno_expansion = [c_rj.get((r, j))*RV.function_stencil_dictionary[WenoConfig.func_points[r*k+j]] for j in range(k)]
            stencil += RV.omega_symbols[r]*sum(eno_expansion)
        RV.reconstructed_expression = stencil
        return


class Weno(Scheme):
    """ Main WENO class. Performs the Jiang-Shu WENO reconstruction procedure. Refer to the reference:
    'Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory schemes for Hyperbolic Conservation
    laws. Shu (1997). The alpha & omega quantities follow the description from this paper."""

    def __init__(self, order, **kwargs):
        """ :arg: int order: Numerical order of the WENO scheme (3,5,...). """
        Scheme.__init__(self, "WenoDerivative", order)
        self.schemetype = "Spatial"
        if order == '5Z':
            # WenoZ only defined for fifth order
            self.k = 3
            self.order = 5
            WT = WenoZ(self.k)
        else:
            self.k = int(0.5*(order+1))
            self.order = order
            WT = WenoJS(self.k)
        JS = JS_smoothness(self.k)
        self.halotype = WenoHalos(self.order)
        self.required_constituent_relations_symbols = {}
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        self.reconstruction_classes = [RightWenoReconstructionVariable('right'), LeftWenoReconstructionVariable('left')]
        # Populate the quantities required by WENO for the left and right reconstruction variable.
        for no, side in enumerate([1, -1]):
            WenoConfig = ConfigureWeno(self.k, side)
            RV = self.reconstruction_classes[no]
            RV.func_points = sorted(set(WenoConfig.func_points))
            RV.stencil_points, RV.function_stencil_dictionary = WenoConfig.generate_symbolic_function_points
            JS.generate_function_smoothness_indicators(RV)
            WT.generate_alphas(RV, WenoConfig)
            WT.generate_omegas(RV, WenoConfig)
            WT.generate_reconstruction(RV, WenoConfig)
            self.reconstruction_classes[no] = RV
        return

    def interpolate_reconstruction_variables(self, derivatives, kernel):
        """ Perform the TENO interpolation on the reconstruction variables.
        arg: list: derivatives: A list of the TENO derivatives to be computed.
        arg: object: kernel: The current computational kernel.
        """
        for d in derivatives:
            pprint(d)
            for rv in d.reconstructions:
                if isinstance(rv, type(self.reconstruction_classes[1])):
                    original_rv = self.reconstruction_classes[1]
                elif isinstance(rv, type(self.reconstruction_classes[0])):
                    original_rv = self.reconstruction_classes[0]
                else:
                    raise ValueError("Reconstruction must be left or right")
                rv.update_quantities(original_rv)
                rv.add_evaluations_to_kernel(kernel)
        return

    def update_constituent_relation_symbols(self, sym, direction):
        """ Function to take the set of required quantities from the constituent relations in symbolic form
         and update the directions in which they are used.
        arg: set: sym: Set of required symbols.
        arg: int: direction: The axis on which WENO is being applied to (x0, x1 ..)."""
        if isinstance(sym, Symbol):
            sym = [sym]
        elif isinstance(sym, list) or isinstance(sym, set):
            pass
        else:
            raise ValueError("The symbol provided should be either a list or symbols")
        for s in sym:

            if isinstance(s, DataSetBase):
                s = s.noblockname
            if s in self.required_constituent_relations_symbols.keys():
                self.required_constituent_relations_symbols[s] += [direction]
            else:
                self.required_constituent_relations_symbols[s] = [direction]
        return

    def generate_constituent_relations_kernels(self, block):
        """ Generates constituent relation kernels on the block.
        arg: object: block: The current block."""
        crs = {}
        for key in self.required_constituent_relations_symbols:
            pprint([key, self.required_constituent_relations_symbols[key]])
            kernel = Kernel(block, computation_name="CR%s" % key)
            kernel.set_grid_range(block)
            for direction in self.required_constituent_relations_symbols[key]:
                kernel.set_halo_range(direction, 0, self.halotype)
                kernel.set_halo_range(direction, 1, self.halotype)
            crs[DataSetBase(str(key))[DataSetBase.location()]] = kernel
        return crs


class EigenSystem(object):
    def __init__(self, eigenvalue, left_ev, right_ev):
        """ Class to hold the routines required by the characteristic decomposition of the Euler equations. The input
        is the eigensystems used to diagonalise the Jacobian in each direction.
        arg: dict: eigenvalue: Dictionary of diagonal matrices containing the eigenvalues in each direction.
        arg: dict: left_ev: Dictionary of matrices containing the left eigenvectors in each direction.
        arg: dict: right_ev: Dictionary of matrices containing the right eigenvectors in each direction. """
        if not isinstance(left_ev, dict) or not isinstance(right_ev, dict) or not isinstance(eigenvalue, dict):
            raise ValueError("Eigen values and eigen vectors should be in dictionary format")
        self.eigen_value = eigenvalue
        self.left_eigen_vector = left_ev
        self.right_eigen_vector = right_ev
        return

    def get_symbols_in_ev(self, direction):
        """ Retrieve the unique symbols present in the eigenvalue matrix for a given direction."""
        return self.eigen_value[direction].atoms(EinsteinTerm).difference(self.eigen_value[direction].atoms(ConstantObject))

    def get_symbols_in_LEV(self, direction):
        """ Retrieve the unique symbols present in the left eigenvector matrix for a given direction."""
        return self.left_eigen_vector[direction].atoms(EinsteinTerm).difference(self.left_eigen_vector[direction].atoms(ConstantObject))

    def get_symbols_in_REV(self, direction):
        """ Retrieve the unique symbols present in the right rigenvector matrix for a given direction."""
        return self.right_eigen_vector[direction].atoms(EinsteinTerm).difference(self.right_eigen_vector[direction].atoms(ConstantObject))

    def generate_grid_variable_ev(self, direction, name):
        """ Create a matrix of eigenvalue GridVariable elements. """
        name = '%s_%d_lambda' % (name, direction)
        return self.symbol_matrix(self.eigen_value[direction], name)

    def generate_grid_variable_REV(self, direction, name):
        """ Create a matrix of right eigenvector GridVariable elements. """
        name = '%s_%d_REV' % (name, direction)
        return self.symbol_matrix(self.right_eigen_vector[direction], name)

    def generate_grid_variable_LEV(self, direction, name):
        """ Create a matrix of left eigenvector GridVariable elements. """
        name = '%s_%d_LEV' % (name, direction)
        return self.symbol_matrix(self.left_eigen_vector[direction], name)

    def symbol_matrix(self, mat, name):
        """ Generic function to populate a matrix of GridVariables for a given name and shape."""
        shape = mat.shape
        symbolic_matrix = zeros(*shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                if mat[i, j]:
                    symbolic_matrix[i, j] = GridVariable('%s_%d%d' % (name, i, j))
        return symbolic_matrix

    def convert_matrix_to_grid_variable(self, mat, name):
        """ Converts the given symbolic matrix to grid variable equivalent.
        arg: Matrix: mat: Symbolic matrix to convert to GridVariable elements.
        arg: str: name: Base name to use for the GridVariables
        returns: mat: The matrix updated to contain GridVariables."""
        syms = list(mat.atoms(EinsteinTerm).difference(mat.atoms(ConstantObject)))
        print type(syms[0])
        print str(syms[0])
        new_syms = [GridVariable('%s_%s' % (name, str(sym))) for sym in syms]
        substitutions = dict(zip(syms, new_syms))
        mat = mat.subs(substitutions)
        return mat

    def generate_equations_from_matrices(self, lhs_matrix, rhs_matrix):
        """ Forms a matrix containing Sympy equations at each element, given two input matrices.
        arg: Matrix: lhs_matrix: The elements of lhs_matrix form the LHS of the output equations.
        arg: Matrix: rhs_matrix: The elements of rhs_matrix form the RHS of the output equations.
        returns: Matrix: equations: A matrix containing an Eq(lhs, rhs) pair in each element."""
        if lhs_matrix.shape != rhs_matrix.shape:
            raise ValueError("Matrices should have the same dimension.")
        equations = zeros(*lhs_matrix.shape)
        for no, v in enumerate(lhs_matrix):
            if rhs_matrix[no] != 0:
                equations[no] = Eq(v, rhs_matrix[no])
        return equations

    def convert_symbolic_to_dataset(self, symbolics, location, direction):
        """ Converts symbolic terms to DataSets.
        arg: object: symbolics: Expression containing Symbols to be converted into DataSets.
        arg: int: location: The integer location to apply to the DataSet.
        arg: int: direction: The integer direction (axis) of the DataSet to apply the location to.
        returns: object: symbolics: The original expression updated to be in terms of DataSets rather than Symbols."""
        symbols = symbolics.atoms(EinsteinTerm).difference(symbolics.atoms(ConstantObject))
        dsets = [DataSetBase(s) for s in symbols]
        loc = DataSetBase.location()
        loc[direction] = loc[direction] + location
        location_datasets = [d[loc] for d in dsets]
        substitutions = dict(zip(symbols, location_datasets))
        return symbolics.subs(substitutions)


class Characteristic(EigenSystem):
    def __init__(self, eigenvalue, left_ev, right_ev):
        """ Class containing the routines required to perform the characteristic decomposition.
        arg: dict: eigenvalue: Dictionary of diagonal matrices containing the eigenvalues in each direction.
        arg: dict: left_ev: Dictionary of matrices containing the left eigenvectors in each direction.
        arg: dict: right_ev: Dictionary of matrices containing the right eigenvectors in each direction. """
        EigenSystem.__init__(self, eigenvalue, left_ev, right_ev)
        self.ev_store = {}
        self.LEV_store = {}
        self.REV_store = {}
        return

    def pre_process(self, direction, derivatives, solution_vector, kernel):
        """ Performs the transformation of the derivatives into characteristic space using the eigensystems provided to Characteristic. Flux splitting is then applied
        to the characteristic variables in preparation for the WENO interpolation. Required quantities are added to pre_process_equations.
        arg: int: direction: Integer direction to apply the characteristic decomposition and WENO (x0, x1, ...).
        arg: list: derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        arg: list: solution_vector: Solution vector from the Euler equations (rho, rhou0, rhou1, rhou2, rhoE) in vector form.
        arg: object: kernel: The current computational kernel."""
        self.direction = direction
        pre_process_equations = []
        required_symbols = self.get_symbols_in_ev(direction).union(self.get_symbols_in_LEV(direction)).union(self.get_symbols_in_REV(direction))
        averaged_suffix_name = 'AVG_%d' % direction
        self.averaged_suffix_name = averaged_suffix_name
        averaged_equations = self.average(required_symbols, direction, averaged_suffix_name)
        pre_process_equations += averaged_equations
        # Eigensystem based on averaged quantities
        avg_LEV_values = self.convert_matrix_to_grid_variable(self.left_eigen_vector[direction], averaged_suffix_name)
        # Grid variables to store averaged eigensystems
        grid_LEV = self.generate_grid_variable_LEV(direction, averaged_suffix_name)
        # Store the grid LEV variables for the characteristic boundary condition
        self.LEV_store[direction] = grid_LEV
        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_LEV, avg_LEV_values))
        # Transform the flux vector and the solution vector to characteristic space
        characteristic_flux_vector, CF_matrix = self.flux_vector_to_characteristic(derivatives, direction, averaged_suffix_name)
        characteristic_solution_vector, CS_matrix = self.solution_vector_to_characteristic(solution_vector, direction, averaged_suffix_name)
        # Can use horner here on characteristic flux
        pre_process_equations += flatten(self.generate_equations_from_matrices(CF_matrix, characteristic_flux_vector))
        pre_process_equations += flatten(self.generate_equations_from_matrices(CS_matrix, characteristic_solution_vector))

        pprint(derivatives)
        # exit()
        for d in derivatives:
            required_symbols = required_symbols.union(d.atoms(DataSetBase))
        self.update_constituent_relation_symbols(required_symbols, direction)

        if hasattr(self, 'flux_split') and self.flux_split:
            print ("Characteristic flux splitting scheme")
            max_wavespeed_matrix, pre_process_equations = self.create_max_characteristic_wave_speed(pre_process_equations, direction)
            # self.split_fluxes(CF_matrix, CS_matrix, max_wavespeed_matrix)
            positive_flux = Rational(1, 2)*(CF_matrix + max_wavespeed_matrix*CS_matrix)
            negative_flux = Rational(1, 2)*(CF_matrix - max_wavespeed_matrix*CS_matrix)
            positive_flux.stencil_points = CF_matrix.stencil_points
            negative_flux.stencil_points = CF_matrix.stencil_points
            self.generate_right_reconstruction_variables(positive_flux, derivatives)
            self.generate_left_reconstruction_variables(negative_flux, derivatives)
        else:
            raise NotImplementedError("Only flux splitting is implemented in characteristic.")
        kernel.add_equation(pre_process_equations)
        return

    def post_process(self, derivatives, kernel):
        """ Transforms the characteristic WENO interpolated fluxes back into real space by multiplying by the right
        eigenvector matrix.
        arg: list: derivatives: The derivatives to perform the characteristic decomposition and WENO on.
        arg: object: kernel: The current computational kernel."""
        post_process_equations = []
        averaged_suffix_name = self.averaged_suffix_name
        avg_REV_values = self.convert_matrix_to_grid_variable(self.right_eigen_vector[self.direction], averaged_suffix_name)
        grid_REV = self.generate_grid_variable_REV(self.direction, averaged_suffix_name)
        self.REV_store[self.direction] = grid_REV
        post_process_equations += flatten(self.generate_equations_from_matrices(grid_REV, avg_REV_values))
        reconstructed_characteristics = Matrix([d.evaluate_reconstruction for d in derivatives])
        reconstructed_flux = grid_REV*reconstructed_characteristics
        reconstructed_work = [d.reconstruction_work for d in derivatives]
        post_process_equations += [Eq(x, y) for x, y in zip(reconstructed_work, reconstructed_flux)]
        kernel.add_equation(post_process_equations)
        return

    def generate_left_reconstruction_variables(self, flux, derivatives):
        if isinstance(flux, Matrix):
            stencil = flux.stencil_points
            for i in range(flux.shape[0]):
                rv = type(self.reconstruction_classes[1])('L_X%d_%d' % (self.direction, i))
                for p in sorted(set(self.reconstruction_classes[1].func_points)):
                    rv.function_stencil_dictionary[p] = flux[i, stencil.index(p)]
                derivatives[i].add_reconstruction_classes([rv])
        return

    def generate_right_reconstruction_variables(self, flux, derivatives):
        if isinstance(flux, Matrix):
            stencil = flux.stencil_points
            for i in range(flux.shape[0]):
                rv = type(self.reconstruction_classes[0])('R_X%d_%d' % (self.direction, i))
                for p in sorted(set(self.reconstruction_classes[0].func_points)):
                    rv.function_stencil_dictionary[p] = flux[i, stencil.index(p)]
                derivatives[i].add_reconstruction_classes([rv])
        return

    def solution_vector_to_characteristic(self, solution_vector, direction, name):
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        solution_vector_stencil = zeros(len(solution_vector), len(stencil_points))
        CS_matrix = zeros(len(solution_vector), len(stencil_points))
        for j, val in enumerate(stencil_points):  # j in fv stencil matrix
            for i, flux in enumerate(solution_vector):
                solution_vector_stencil[i, j] = incr_dset(flux, direction, val)
                CS_matrix[i, j] = GridVariable('CS_%d%d' % (i, j))
        grid_LEV = self.generate_grid_variable_LEV(direction, name)
        characteristic_solution_stencil = grid_LEV*solution_vector_stencil
        characteristic_solution_stencil.stencil_points = stencil_points
        CS_matrix.stencil_points = stencil_points
        return characteristic_solution_stencil, CS_matrix

    def flux_vector_to_characteristic(self, derivatives, direction, name):
        fv = []
        for d in derivatives:
            if d.get_direction[0] != direction:
                raise ValueError("Derivatives provided for flux vector are not homogeneous in direction.")
            fv += [d.args[0]]
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        flux_stencil = zeros(len(fv), len(stencil_points))
        CF_matrix = zeros(len(fv), len(stencil_points))
        for j, val in enumerate(stencil_points):  # j in fv stencil matrix
            for i, flux in enumerate(fv):
                flux_stencil[i, j] = incr_dset(flux, direction, val)
                CF_matrix[i, j] = GridVariable('CF_%d%d' % (i, j))
        grid_LEV = self.generate_grid_variable_LEV(direction, name)
        characteristic_flux_stencil = grid_LEV*flux_stencil
        characteristic_flux_stencil.stencil_points = stencil_points
        CF_matrix.stencil_points = stencil_points
        return characteristic_flux_stencil, CF_matrix


class LLFCharacteristic(Characteristic):
    """ This class contains the Local Lax-Fedrich scheme
    """

    def __init__(self, eigenvalue, left_ev, right_ev, order, ndim, averaging=None):
        Characteristic.__init__(self, eigenvalue, left_ev, right_ev)
        # Weno.__init__(self, order)
        if averaging is None:
            self.average = SimpleAverage([0, 1]).average
        else:
            self.average = averaging.average
        self.ndim = ndim
        self.flux_split = True
        return

    def create_max_characteristic_wave_speed(self, pre_process_equations, direction):
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        ev = self.eigen_value[direction]
        out = zeros(*ev.shape)
        for p in stencil_points:
            location_ev = self.convert_symbolic_to_dataset(ev, p, direction)
            for no, val in enumerate(location_ev):
                out[no] = Max(out[no], Abs(val))
        max_wave_speed = self.generate_grid_variable_ev(direction, 'max')
        ev_equations = self.generate_equations_from_matrices(max_wave_speed, out)
        ev_equations = [x for x in ev_equations if x != 0]
        pre_process_equations += ev_equations
        max_wave_speed = diag(*([max_wave_speed[i, i] for i in range(ev.shape[0])]))
        return max_wave_speed, pre_process_equations

    def discretise(self, type_of_eq, block):
        """ This is the place where the logic of vector form of equations are implemented.
        Find physical fluxes by grouping derivatives by direction --> in central, copy over
        Then the physical fluxes are transformed to characteristic space ---> a function in Characteristic
        Find max lambda and update f+ and f- ------> This would be in this class (GLF)
        For each f+ and f-, find f_hat of i+1/2, i-1/2, (L+R) are evaluated  ----> Function in WENO scheme, called from in here
        flux at i+1/2 evaluated -- > Function in WENO scheme
        Then WENO derivative class is instantiated with the flux at i+1/2 array --> Function in WENO scheme, called from in here
        Final derivatives are evaluated from Weno derivative class --> Using WD.discretise
        """
        if isinstance(type_of_eq, SimulationEquations):
            # residual_kernel = Kernel(block, computation_name = "WenoResidual")
            eqs = flatten(type_of_eq.equations)
            grouped = self.group_by_direction(eqs)
            all_derivatives_evaluated_locally = []

            reconstruction_halos = self.reconstruction_halotype(self.order, reconstruction=True)
            for key, derivatives in grouped.iteritems():
                all_derivatives_evaluated_locally += derivatives
                for no, deriv in enumerate(derivatives):
                    deriv.create_reconstruction_work_array(block)
                kernel = Kernel(block, computation_name="%s_reconstruction_%d_direction" % (self.__class__.__name__, key))
                kernel.set_grid_range(block)
                # WENO reconstruction should be evaluated for extra point on each side
                kernel.set_halo_range(key, 0, reconstruction_halos)
                kernel.set_halo_range(key, 1, reconstruction_halos)
                self.pre_process(key, derivatives, flatten(type_of_eq.time_advance_arrays), kernel)
                self.interpolate_reconstruction_variables(derivatives, kernel)
                block.set_block_boundary_halos(key, 0, self.halotype)
                block.set_block_boundary_halos(key, 1, self.halotype)
                type_of_eq.Kernels += [kernel]
            if grouped:
                type_of_eq.Kernels += [self.evaluate_residuals(block, eqs, all_derivatives_evaluated_locally)]
            constituent_relations = self.generate_constituent_relations_kernels(block)
            return constituent_relations

    def evaluate_residuals(self, block, eqns, local_ders):
        residue_eq = []
        for eq in eqns:
            substitutions = {}
            for d in eq.rhs.atoms(Function):
                if d in local_ders:
                    substitutions[d] = d._discretise_derivative(block)
                else:
                    substitutions[d] = 0
            residue_eq += [Eq(eq.residual, eq.residual + eq.rhs.subs(substitutions))]
        residue_kernel = Kernel(block, computation_name="%s Residual" % self.__class__.__name__)
        residue_kernel.set_grid_range(block)
        residue_kernel.add_equation(residue_eq)
        return residue_kernel


class SimpleAverage(object):
    def __init__(self, locations):
        self.locations = locations
        return

    def average(self, functions, direction, name_suffix):
        """Performs simple average.
        arg: functions: List of function (Symbols) to apply averaging on.
        arg: locations: Relative index used for averaging (e.g. [0,1] for i and i+1)
        arg: direction: Axis of the dataset on which the location should be applied.
        arg; name_suffix: Name to be appended to the functions. """
        avg_equations = []
        for f in functions:
            name = f.get_base()
            d = DataSetBase(name)
            loc = [0 for i in range(DataSetBase.block.ndim)]
            loc1 = loc[:]
            loc2 = loc[:]
            loc1[direction] = loc[direction] + self.locations[0]
            loc2[direction] = loc[direction] + self.locations[1]
            a = d[loc1]
            b = d[loc2]
            avg_equations += [Eq(GridVariable('%s_%s' % (name_suffix, name)), (a+b)/2)]
        return avg_equations


class RoeAverage(object):
    def __init__(self, locations, ndim, physics=None):
        self.locations = locations
        self.ndim = ndim
        if not physics:
            self.NS = NSphysics(ndim)
        else:
            self.NS = physics
        return

    def get_dsets(self, dset_base):
        loc = dset_base.location()
        loc1, loc2 = loc[:], loc[:]
        loc1[self.direction], loc2[self.direction] = loc[self.direction]+self.locations[0], loc[self.direction]+self.locations[1]
        dset1, dset2 = dset_base[loc1], dset_base[loc2]
        return dset1, dset2

    def average(self, functions, direction, name_suffix):
        self.direction = direction
        evaluations = []
        names = []
        # Averaged density rho_hat = sqrt(rho_L*rho_R)
        rho_base = self.NS.density().base
        rho_L, rho_R = self.get_dsets(rho_base)
        evaluations += [sqrt(rho_L*rho_R)]
        grid_vars = [GridVariable('%s_%s' % (name_suffix, rho_base.simplelabel()))]
        # Store inverse factor 1/(sqrt(rho_R)+sqrt(rho_L))
        grid_vars += [GridVariable('%s_%s' % (name_suffix, 'inv_rho'))]
        evaluations += [(sqrt(rho_R)+sqrt(rho_L))**(-1)]

        # Average velocity components
        for i, velocity in enumerate(self.NS.velocity()):
            velocity_base = velocity.base
            vel_L, vel_R = self.get_dsets(velocity_base)
            names += [velocity_base.simplelabel()]
            evaluations += [(sqrt(rho_L)*vel_L+sqrt(rho_R)*vel_R)*grid_vars[1]]
            grid_vars += [GridVariable('%s_%s' % (name_suffix, velocity_base.simplelabel()))]
        # Averaged kinetic energy in terms of u0, u1, u2 grid variables
        KE = Rational(1, 2)*sum([u**2 for u in grid_vars[2:]])

        # Calcualte enthalpy h = rhoE + P/rho
        pressure_base = self.NS.pressure().base
        P_L, P_R = self.get_dsets(pressure_base)
        rhoE_base = self.NS.total_energy().base
        rhoE_L, rhoE_R = self.get_dsets(rhoE_base)
        H_L = rhoE_L + P_L/rho_L
        H_R = rhoE_R + P_R/rho_R
        roe_enthalpy = (sqrt(rho_L)*H_L+sqrt(rho_R)*H_R)*grid_vars[1]
        # Average speed of sound
        a_base = self.NS.speed_of_sound().base
        roe_a = sqrt((self.NS.specific_heat_ratio()-1)*(roe_enthalpy - KE))
        evaluations += [roe_a]
        grid_vars += [GridVariable('%s_%s' % (name_suffix, a_base.simplelabel()))]

        return [Eq(x, y) for (x, y) in zip(grid_vars, evaluations)]


class LLFWeno(LLFCharacteristic, Weno):
    def __init__(self, eigenvalue, left_ev, right_ev, order, ndim, averaging=None):
        LLFCharacteristic.__init__(self, eigenvalue, left_ev, right_ev, order, ndim, averaging)
        Weno.__init__(self, order)
        return

    def reconstruction_halotype(self, order, reconstruction=True):
        return WenoHalos(order, reconstruction)

    def group_by_direction(self, eqs):
        """ Groups the input equations by the direction (x0, x1, ...) they depend upon.
        arg: list: eqs: List of equations to group by direction.
        returns: dict: grouped: Dictionary of {direction: equations} key, value pairs for equations grouped by direction."""
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(WenoDerivative))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        # TODO: check for size of grouped items
        return grouped

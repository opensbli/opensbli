from sympy import *
from .opensblifunctions import *
from opensbli.core import *
from .grid import GridVariable

class Scheme(object):
    """ A numerical discretisation scheme. """
    def __init__(self, name, order):
        """ Initialise the scheme.

        :arg str name: The name of the scheme.
        :arg int order: The order of the scheme.
        """
        self.name = name
        self.order = order
        return

class WenoHalos(object):
    def __init__(self, order):
        # Check for the boundary types in the blocks and set the halo points
        #self.halos = [[-scheme.order, scheme.order] for dim in range(block.ndim)]
        k = int(0.5*(order+1))
        # self.halos = [(-k, k+1)for dim in range(ndim)]
        self.halos = [-k, k+1]
        return
    def get_halos(self, side):
        return self.halos[side]

class ConfigureWeno(object):
    def __init__(self, k, side):
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
        
        return
    @property
    def generate_symbolic_function_points(self):
        ndim = DataSet.dimensions
        DataSet.dimensions = 1
        f = DataSet('f')
        symbolic_functions = []
        symbolic_function_dictionary = {}
        for p in set(self.func_points):
            symbolic_functions.append(f.get_location_dataset([p]))
            symbolic_function_dictionary[p] = symbolic_functions[-1]
        DataSet.dimensions = ndim
        return symbolic_functions, symbolic_function_dictionary

    def generate_eno_coefficients(self, k):
        """ Generates the c_rj ENO reconstruction coefficients given in Table 2.1
        of 'Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes
        for Hyperbolic Conservation Laws' by Shu(1997). Computation is of equation (2.21).

        :arg int k: WENO coefficient k, equal to scheme order = 2k - 1.
        :arg int side: Reconstruction side. -1 for left, 1 for right.
        :returns: dict c_rj: Dictionary in the form (r,j) : ENO coefficient c_rj.
        """
        side = self.side
        if side == 1:
            d = 1
        elif side == -1:
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
                        for q in [x for x in range(k+1) if (x != m and x!= l)]:
                            top_product = r - q + d
                            top_product_total *= top_product
                        bottom_product *= m - l
                        top_sum += top_product_total
                    c_rj_sum += Rational(top_sum, bottom_product)
                c_rj[(r,j)] = c_rj_sum
        return c_rj

    def generate_optimal_weights(self):
      """
      """
      # Store reconstruction coefficients to form the sum
      k, c_rj, c2_rj = self.k, self.c_rj, self.c2_rj
      opt_weights = [Symbol('d_%d' % r) for r in range(k)]
      r_values = [[k-1-j for j in range(m+1)][::-1] for m in range(k)]
      equations = []
      for i, indices in enumerate(r_values):
        variables = [opt_weights[r] for r in indices]
        coeffs = [c_rj[(r,j)] for j, r in enumerate(indices)]
        rhs = c2_rj[(k-1, i)]
        v = [variables[i]*coeffs[i] for i in range(len(indices))]
        equations += [sum(v) - rhs]
      # Solve the linear system to get the coefficients
      solution = solve(equations, opt_weights)
      coeff_dict = {}
      for i in range(k):
        coeff_dict[(0,i)] = solution[opt_weights[i]]
      return coeff_dict

    def generate_left_right_points(self):
        """ Populate the function evaluation points for left and right reconstruction.
        args: None
        returns: list: fn_points
        """
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
        self.k = k
        self.smooth_coeffs = self.generate_smoothness_coefficients()
        return

    def generate_smoothness_coefficients(self):
      """Extracts the JS smoothness coefficients."""
      k = self.k
      smooth_coeffs = {}
      x, dx = Symbol('x'), Symbol('dx')
      for r in range(k):
        # Generate x, y values to create interpolating polynomial
        dx_values = [dx*i for i in range(-r, k+1-r)]
        nodes = [Symbol('fn[i%+d]' % i) for i in range(-r,k-r)]
        funcs = [0]
        for i in range(k):
          funcs.append(funcs[-1] + dx*nodes[i])
        # Perform Lagrange interpolation
        lagrange_interp = interpolating_poly(k+1, x, dx_values, funcs).diff(x)
        # Loop over the derivatives of the interpolating polynomial
        total = 0
        for l in range(1,k):
          q = (lagrange_interp.diff(x, l))**2
          # Perform the integration and multiply by h^(2*l-1) over cell
          q = integrate(q.as_poly(x), x) * dx**(2*l-1)
          total += (q.subs(x, dx) - q.subs(x,0))
        done = []
        # Save the coefficients of the smoothness indicator
        for m in range(0, 2*k-2):
          for n in range(0, 2*k-2):
            func_product = Symbol('fn[i%+d]' % (-r+m))*Symbol('fn[i%+d]' % (-r+n))
            if func_product not in done:
              c = total.coeff(func_product)
              if c != 0:
                smooth_coeffs[(r,m,n)] = c
              done.append(func_product)
      return smooth_coeffs

    def generate_function_smoothness_indicators(self, fn):
        if isinstance(fn, LeftReconstructionVariable):
            shift = 1
        elif isinstance(fn, RightReconstructionVariable):
            shift = 0
        k = self.k
        # Compute the smoothness indicator and alpha
        for r in range(k):
            fn.smoothness_symbols += [Symbol('beta_%d' % r)]
            local_smoothness = 0
            for m in range(0, 2*k-2): # WARNING: change to (0, 2*k -1) if this is now broken.
                for n in range(0, 2*k-2):
                    beta = self.smooth_coeffs.get((r,m,n))
                    if beta != None:
                        shift_indices = [-r+m+shift, -r+n+shift]
                        func_product = fn.function_stencil_dictionary[shift_indices[0]]*fn.function_stencil_dictionary[shift_indices[1]]
                        local_smoothness += beta*func_product
            fn.smoothness_indicators += [local_smoothness]
        return

class ReconstructionVariable(object):
    def __init__(self, name):
        self.smoothness_indicators = []
        self.smoothness_symbols = []
        self.alpha_evaluated = []
        self.alpha_symbols = []
        self.omega_evaluated = []
        self.omega_symbols = []
        return
    @property
    def create_equations(self):
        all_symbols = self.smoothness_symbols + self.alpha_symbols + self.omega_symbols
        all_evaluations = self.smoothness_indicators + self.alpha_evaluated + self.omega_evaluated

        all_equations = (zip(all_symbols, all_evaluations))
        return all_equations

    def create_function_equations(self, drv):
        all_eqs = drv.create_equations
        old_eq = all_eqs[:]
        subs_dict = {}
        for p in drv.func_points:
            subs_dict[drv.function_stencil_dictionary[p]] = self.function_stencil_dictionary[p]
        pprint(subs_dict)
        for no, eq in enumerate(all_eqs):
            all_eqs[no] = tuple([eq[0].subs(subs_dict), eq[1].subs(subs_dict)])
        for no, eq in enumerate(all_eqs):
            pprint(Eq(eq[0], eq[1]))
            pprint(Eq(old_eq[no][0], old_eq[no][1]))
        # pprint(all_eqs)
        return

class LeftReconstructionVariable(ReconstructionVariable):
    def __init__(self, name):
        ReconstructionVariable.__init__(self, name)
        return

class RightReconstructionVariable(ReconstructionVariable):
    def __init__(self, name):
        ReconstructionVariable.__init__(self, name)
        return


class Weno(Scheme):
    """ Main WENO class."""
    def __init__(self, order, **kwargs):
        """
        :returns: None
        """
        Scheme.__init__(self, "WenoDerivative", order)
        self.eps = Symbol('epsilon')
        k = int(0.5*(order+1))
        JS = JS_smoothness(k)
        self.schemetype = "Spatial"
        self.k = k
        self.halotype = WenoHalos(order)
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        smooth_coeffs = JS.smooth_coeffs
        reconstruction_classes = [LeftReconstructionVariable('left'), RightReconstructionVariable('right')]
        for no, side in enumerate([-1, 1]):

            WenoConfig = ConfigureWeno(self.k, side)
            rc = reconstruction_classes[no]
            rc.func_points = sorted(set(WenoConfig.func_points))
            # pprint(WenoConfig.func_points)

            rc.stencil_points, rc.function_stencil_dictionary = WenoConfig.generate_symbolic_function_points
            # rc.create_stencil_dictionary
            JS.generate_function_smoothness_indicators(rc)
            
            self.generate_alphas(rc, WenoConfig)
            self.generate_omegas(rc, WenoConfig)
            self.generate_reconstruction(rc, WenoConfig)
            reconstruction_classes[no] = rc

        a = LeftReconstructionVariable('test')
        a.function_stencil_dictionary = {}
        t = DataSet('e')
        fn = (DataSet('p')+DataSet('rhoE'))*DataSet('u0')
        fn = DataSet('rhoE')*DataSet('u0')
        for k in reconstruction_classes[0].func_points:
            loc = t.location
            loc[0] = loc[0] + k
            local_dict = {}
            for d in fn.atoms(DataSet):
                local_dict[d] = d.get_location_dataset(loc)
            a.function_stencil_dictionary[k] = fn.subs(local_dict)
        pprint(a.function_stencil_dictionary)
        reconstruction_classes[0].create_equations
        a.create_function_equations(reconstruction_classes[0])
        exit()
        return

    def generate_weno_reconstruction(self, fn):

        return

    def generate_alphas(self, fn, WenoConfig):
        for r in range(self.k):
            fn.alpha_symbols += [Symbol('alpha_%d' % r)]
            fn.alpha_evaluated.append(WenoConfig.opt_weights.get((0,r))/(self.eps+fn.smoothness_symbols[r])**2)
        return

    def generate_omegas(self, fn, WenoConfig):
        for r in range(self.k):
            fn.omega_symbols += [Symbol('omega_%d' % r)]
            fn.omega_evaluated += [fn.alpha_symbols[r]/sum(fn.alpha_symbols)]
        return

    def generate_reconstruction(self, fn, WenoConfig):
        stencil = 0
        k, c_rj = self.k, WenoConfig.c_rj
        for r in range(k):
            eno_expansion = [c_rj.get((r,j))*fn.function_stencil_dictionary[WenoConfig.func_points[r*k+j]] for j in range(k)]
            stencil += fn.omega_symbols[r]*sum(eno_expansion)

        fn.reconstruction = stencil
        return


class EigenSystem(object):
    def __init__(self, eigenvalue, left_ev, right_ev):
        if not isinstance(left_ev, dict) or not isinstance(right_ev, dict) or not isinstance(eigenvalue, dict):
            raise ValueError("Eigen values and eigen vectors should be in dictionary format")
        self.eigen_value = eigenvalue
        self.left_eigen_vector = left_ev
        self.right_eigen_vector = right_ev
        return

    def get_symbols_in_ev(self, direction):
        return self.eigen_value[direction].atoms(Symbol)

    def get_symbols_in_LEV(self, direction):
        return self.left_eigen_vector[direction].atoms(Symbol)

    def get_symbols_in_REV(self, direction):
        return self.right_eigen_vector[direction].atoms(Symbol)

    def generate_grid_variable_ev(self, direction, name):
        name = '%s_%d_lambda' % (name, direction)
        return self.symbol_matrix(self.eigen_value[direction].shape, name)

    def generate_grid_variable_REV(self, direction, name):
        name = '%s_%d_REV' % (name, direction)
        return self.symbol_matrix(self.right_eigen_vector[direction].shape, name)

    def generate_grid_variable_LEV(self, direction, name):
        name = '%s_%d_LEV' % (name, direction)
        return self.symbol_matrix(self.left_eigen_vector[direction].shape, name)

    def symbol_matrix(self, shape, name):
        symbolic_matrix = zeros(*shape)
        for i in range(shape[0]):
            for j in range(shape[1]):
                symbolic_matrix[i,j] = GridVariable('%s_%d%d'%(name,i,j))
        return symbolic_matrix

    def convert_matrix_to_grid_variable(self, mat, name):
        """Converts the given matrix function to grid variable equivalent
        """
        syms = list(mat.atoms(Symbol))
        new_syms = [GridVariable('%s_%s'%(name,str(sym))) for sym in syms]
        substitutions = dict(zip(syms, new_syms))
        mat = mat.subs(substitutions)
        return mat

    def generate_equations_from_matrices(self, m1, m2):

        if m1.shape != m2.shape:
            raise ValueError("Matrices should have the same dimension.")
        equations = zeros(*m1.shape)
        for no, v in enumerate(m1):
            if m2[no] != 0:
                equations[no] = Eq(v, m2[no])
        return equations

    def increment_dataset(self, expression, direction, value):
        """
        Increments a dataset by the given increment in the direction passed to the characteristic routine.
        """
        for dset in expression.atoms(DataSet):
            loc = dset.location[:]
            loc[direction] = loc[direction] + value
            new_dset =  dset.get_location_dataset(loc)
            expression = expression.replace(dset, new_dset)
        return expression





class Characteristic(EigenSystem):
    def __init__(self,  eigenvalue, left_ev, right_ev):
        """ This holds the logic for the reconstruction procedure"""
        self.is_vector_type = True
        self.leftRight = [True, True]
        EigenSystem.__init__(self, eigenvalue, left_ev, right_ev)
        return

    def group_by_direction(self, eqs):
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(WenoDerivative))
        all_WDS = list(set(all_WDS))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped


    def pre_process(self, direction, derivatives, solution_vector):
        """
        """
        pre_process_equations = []
        required_symbols = self.get_symbols_in_ev(direction).union(self.get_symbols_in_LEV(direction)).union(self.get_symbols_in_REV(direction))
        averaged_suffix_name = 'AVG'
        averaged_equations = self.average(required_symbols, direction, averaged_suffix_name)
        # Eigensystem based on averaged quantities
        avg_EV_values = self.convert_matrix_to_grid_variable(self.eigen_value[direction], averaged_suffix_name)
        avg_REV_values = self.convert_matrix_to_grid_variable(self.right_eigen_vector[direction], averaged_suffix_name)
        avg_LEV_values = self.convert_matrix_to_grid_variable(self.left_eigen_vector[direction], averaged_suffix_name)
        # Grid variables to store averaged eigensystems
        grid_EV = self.generate_grid_variable_ev(direction, averaged_suffix_name)
        grid_LEV = self.generate_grid_variable_LEV(direction, averaged_suffix_name)
        grid_REV = self.generate_grid_variable_REV(direction, averaged_suffix_name)

        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_EV ,avg_EV_values))
        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_LEV ,avg_LEV_values))
        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_REV ,avg_REV_values))
        pre_process_equations = averaged_equations + [x for x in pre_process_equations if x != 0]
        # Transform the flux vector and the solution vector to characteristic space
        characteristic_flux_vector = self.flux_vector_to_characteristic(derivatives)
        characteristic_solution_vector = self.solution_vector_to_characteristic(solution_vector, direction)
        pprint(characteristic_flux_vector[:,0])
        print "-----"
        pprint(characteristic_solution_vector[0,0])
        # Here need to equate the symbolic CF, CS, alpha with the entries in characteristic flux/solution vectors to form the flux splitting
        CF_matrix, CS_matrix, CW_speed_matrix = self.create_characteristic_symbol_matrices()
        CF_matrix.stencil_points = characteristic_flux_vector.stencil_points
        pprint(CF_matrix.stencil_points)
        pprint(CF_matrix)
        pprint(CS_matrix)
        exit()

        #######
        side = -1
        index = 0
        left_WENO_reconstruction = self.generate_weno(direction, side, index)
        side = 1
        right_WENO_reconstruction = self.generate_weno(direction, side, index)
        exit()

        if hasattr(self, 'flux_split') and self.flux_split:
            print ("Characteristic flux splitting scheme")
            # If flux splitting scheme then evaluate characteristic solution
            interpolated_solution = self.apply_weno(char_solution_evaluation)
            pprint(char_solution_evaluation)
            exit()
            # Evaluate split parameter, this can be an equation or a Kernel
            # Perform splitting of fluxes
            # Do the interpolation
        else:
            print ("Characteristic WENO solution without splitting fluxes")
            pass

        """ TODO
        Add flux vector characteristic n solution vector characteristic evaluations to pre_process_equations
        Find the maximum lambda for the stencil
        Create left flux and right flux (f(u) + alpha *u and f(u) - alpha*u
        perform weno on these
        """
        exit()

        return

    def create_characteristic_symbol_matrices(self):
        CF_matrix = zeros(self.ndim+2, self.order+1)
        CS_matrix = zeros(self.ndim+2, self.order+1)
        CW_speed_matrix = zeros(self.ndim+2, self.ndim+2)
        for i in range(self.ndim+2):
            for j in range(self.order+1):
                CF_matrix[i,j] = Symbol('CF_%d%d' % (i,j))
                CS_matrix[i,j] = Symbol('CS_%d%d' %(i,j))
            CW_speed_matrix[i,i] = Symbol('WS_%i%i' % (i,i))
        return CF_matrix, CS_matrix, CW_speed_matrix


    def apply_weno(self, char_solution_evaluation):

        return

    def solution_vector_to_characteristic(self, solution_vector, direction):
        stencil_points = sorted(list(set(self.WenoConfigL.func_points + self.WenoConfigR.func_points)))
        print "stencil points are: ", stencil_points
        solution_vector_stencil = zeros(len(solution_vector), len(stencil_points))
        for j, val in enumerate(stencil_points): # j in fv stencil matrix
            for i, flux in enumerate(solution_vector):
                solution_vector_stencil[i,j] = self.increment_dataset(flux, direction, val)
        grid_LEV = self.generate_grid_variable_LEV(direction, 'AVG')
        characteristic_solution_stencil = grid_LEV*solution_vector_stencil
        characteristic_solution_stencil.stencil_points = stencil_points
        return characteristic_solution_stencil
    def flux_vector_to_characteristic(self, derivatives):
        directions = []
        fv = []
        for d in derivatives:
            directions += d.get_direction
            fv += [d.args[0]]
        if len(set(directions)) != 1:
            raise ValueError("Derivatives provided for flux vector are not homogeneous in direction.")
        direction = directions[0]
        stencil_points = sorted(list(set(self.WenoConfigL.func_points + self.WenoConfigR.func_points)))
        flux_stencil = zeros(len(fv), len(stencil_points))
        for j, val in enumerate(stencil_points): # j in fv stencil matrix
            for i, flux in enumerate(fv):
                flux_stencil[i,j] = self.increment_dataset(flux, direction, val)
        grid_LEV = self.generate_grid_variable_LEV(direction, 'AVG')
        characteristic_flux_stencil = grid_LEV*flux_stencil
        characteristic_flux_stencil.stencil_points = stencil_points
        return characteristic_flux_stencil

    def create_dictionary_interpolations(self, interpolated):
        dictionary_interp = {}
        direction = set()
        for i in interpolated:
            dictionary_interp[i.variable] = i
            direction.add(i.direction)
        if len(direction) > 1:
            raise ValueError("Something is wrong")
        direction = list(direction)[0]
        return dictionary_interp, direction

    def post_process(self, interpolated):
        """ The left and right interpolations of the characteristic variables are returned back.
        The process is first evaluate the averaged state, left and right eigen vectors and eigen
        values
        After that write the equations for the interpolations
        """
        self.post_process_equations = []
        direction = self.direction
        # Equations for the evaluation of interpolations and store the reconstructions
        Interpolated_left_characteristic_flux = [] # Should be a list as we do not want to change the order
        Interpolated_right_characteristic_flux = []
        Interpolated_right_solution = []
        Interpolated_left_solution = []
        dictionary, key = self.create_dictionary_interpolations(interpolated)

        temp_dictionary = {} # Temporary dictionary to store the reconstructed_symbols
        # Do the evaluations for the interpolated quantities
        for val in interpolated:
            self.post_process_equations, reconstructed_symbols = val.evaluate_interpolation(self.post_process_equations)
            temp_dictionary[val.variable] = reconstructed_symbols

        # The new naming uplus will be minus
        self.right_interpolated = []
        self.left_interpolated = []
        # Identify the right/left reconstructions
        for val in self.u_plus:
            self.right_interpolated += temp_dictionary[val][1]
        for val in self.u_minus:
            self.left_interpolated += temp_dictionary[val][-1]

        # Construct the final flux using the reconstruction placeholders
        # and transform back to the physical space using the symbolic right eigenvector matrix REV.
        final_flux = self.REV_symbolic*(Matrix(self.right_interpolated)) + \
            self.REV_symbolic*(Matrix(self.left_interpolated))
        final_equations = self.pre_process_equations + self.post_process_equations
        # pprint(self.pre_process_equations)
        # print "##############################"
        # pprint(self.post_process_equations)
        # from .latex import *

        # fname = './hyper_eq.tex'
        # latex = LatexWriter()
        # latex.open('./hyper_eq.tex')
        # metadata = {"title": "Hyperbolic Flux equations", "author": "David", "institution": ""}
        # latex.write_header(metadata)
        # s = "Pre_process equations: "
        # latex.write_string(s)
        # for eq in self.pre_process_equations:
        #     latex.write_expression(eq)
        # s = "Post_process equations: "
        # latex.write_string(s)
        # for eq in self.post_process_equations:
        #     latex.write_expression(eq)
        # latex.write_footer()
        # latex.close()
        return final_equations, final_flux

class LLFCharacteristic(Characteristic, Weno):
    """ This class contains the Local Lax-Fedrich scheme
    """
    def __init__(self,  eigenvalue, left_ev, right_ev, order, ndim, averaging=None):
        Characteristic.__init__(self, eigenvalue, left_ev, right_ev)
        Weno.__init__(self, order)
        if averaging == None:
            self.average = SimpleAverage([0,1]).average
        else:
            self.average = averaging.average
        self.ndim = ndim
        self.flux_split = True
        return

    def split_fluxes(self, char_flux, char_solution):

        return

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
        #print "Here "
        from .opensbliequations import *
        if isinstance(type_of_eq, SimulationEquations):
            eqs = flatten(type_of_eq.equations)
            grouped = self.group_by_direction(eqs)
            for key, value in grouped.iteritems():
                print(key, value)
                for v in value:
                    v.update_work(block)
                    print v.__dict__, v
                self.pre_process(key, value, flatten(type_of_eq.time_advance_arrays))
            # pprint(type_of_eq.__dict__)

            exit()
        return

    def group_by_direction(self, eqs):
        all_WDS = []
        grouped = {}
        for eq in eqs:
            local_wds = list(eq.atoms(WenoDerivative))
            all_WDS = list(eq.atoms(WenoDerivative))
            for wd in all_WDS:
                direction = wd.get_direction[0]
                if direction in grouped.keys():
                    grouped[direction] += [wd]
                else:
                    grouped[direction] = [wd]
        return grouped

    def interp_functions(self, key):
        local_equations = []
        # Symbols for absolute max eigenvalues
        f = lambda x:x.subs(x, abs(x))
        centre = self.eigenvalue_names[0].applyfunc(f)
        diagonal_entries = [e_val for e_val in centre if e_val != 0]
        max_lambda_names = [GridVariable('max_lambda')]
        local_equations = [Eq(max_lambda_names[0], Max(*(diagonal_entries)))]
        lambda_max_matrix = diag(*([max_lambda_names[0]]*centre.shape[0]))

        flux_vector = Matrix(self.vector_notation[key])
        time_vector = Matrix(self.vector_notation[CoordinateObject('t')])
        ### Here multiply the LEV symbolic matrix (LEV00 etc) by the flux and solution vectors
        ### Transformation to characteristic space is done here THIS IS PART (c) of the WENO procedure 2.10
        flux_conserve = self.LEV_symbolic*flux_vector
        time_vector = self.LEV_symbolic*Matrix(self.vector_notation[CoordinateObject('t')])
        # Perform the flux splitting by finding max absolute eigenvalues
        positive = []
        negative = []
        max_lam_vec = []
        for i in range(len(flux_vector)):
            lm2EV = self.eigenvalue_names[-2][i,i]
            l_EV = self.eigenvalue_names[-1][i,i]
            c_EV = self.eigenvalue_names[0][i,i]
            r_EV = self.eigenvalue_names[1][i,i]
            r2EV = self.eigenvalue_names[2][i,i]
            r3EV = self.eigenvalue_names[3][i,i] # 3rd not used?
            # The solution for if local lambda == 0
            max_lambda = Max(Abs(c_EV),Abs(r_EV), Abs(l_EV), Abs(r2EV),Abs(lm2EV))
            max_lam_vec += [GridVariable('max_lambda_%d'%i)]
            # Equates the above grid variable max_lambda_0 to the max_lambda equation expression for this component of flux vector
            # Adds to the massive pre_process list
            self.pre_process_equations += [Eq(max_lam_vec[-1],max_lambda )]
            # This part does f+ = 0.5*(u(i) + MAX(alpha?)*f(u(i)))
            # and f- = 0.5*(u(i) - MAX(alpha?)*f(u(i)))
            positive += [Rational(1,2)*(flux_conserve[i] + max_lam_vec[-1]*time_vector[i])]
            negative += [Rational(1,2)*(flux_conserve[i] - max_lam_vec[-1]*time_vector[i])]

        self.max_lam_vec = max_lam_vec

        # The fluxes in characteristic form are stored to self.u_plus/minus
        self.u_plus = positive
        self.u_minus = negative
        ## At this point the FULL flux equations in both u_plus/minus are passed to WENO to do the interpolations, with their key

        ## weno procedure is applied in update_WenoSolutionType, which is called on the pre_processed equations in the main calling computations part
        ## This is done after WenoSolutionType is called updating which flux terms need WENO and the side to apply WENO to
        ## All of the below is setting up what weno procedures are needed, the actual weno is applied by update_WenoSolutionType.
        required_interpolations = []
        leftRight = [False, True]
        for val in self.u_plus:
            temp = WenoSolutionType(val,leftRight)
            temp.direction = key
            temp.direction_index = self.direction_index
            required_interpolations += [temp]
        leftRight = [True, False]
        for val in self.u_minus:
            temp = WenoSolutionType(val,leftRight)
            temp.direction = key
            temp.direction_index = self.direction_index
            required_interpolations += [temp]
        return required_interpolations


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
            d = DataSet(str(f))
            loc = d.location[:]
            loc1 = loc[:]
            loc2 = loc[:]
            loc1[direction] = loc[direction] + self.locations[0]
            loc2[direction] = loc[direction] + self.locations[1]
            a = d.get_location_dataset(loc1)
            b = d.get_location_dataset(loc2)
            avg_equations += [Eq(GridVariable('%s_%s' % (name_suffix, str(f))), (a+b)/2)]
        return avg_equations
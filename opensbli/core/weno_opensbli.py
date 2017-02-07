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
    @property
    def create_equations(self):
        all_symbols = self.smoothness_symbols + self.alpha_symbols + self.omega_symbols
        all_evaluations = self.smoothness_indicators + self.alpha_evaluated + self.omega_evaluated

        all_equations = (zip(all_symbols, all_evaluations))
        pprint(all_equations)
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

    def update(self, original):
        print self.name
        self.smoothness_symbols += [GridVariable('%s%s' % (s, self.name)) for s in original.smoothness_symbols]
        self.alpha_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.alpha_symbols]
        self.omega_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.omega_symbols]

        subs_dict = dict(zip(original.smoothness_symbols+original.alpha_symbols+original.omega_symbols, self.smoothness_symbols+self.alpha_symbols+self.omega_symbols))
        fn_subs_dict = {}
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
        self.reconstruction_classes = [LeftReconstructionVariable('left'), RightReconstructionVariable('right')]
        for no, side in enumerate([-1, 1]):

            WenoConfig = ConfigureWeno(self.k, side)
            rc = self.reconstruction_classes[no]
            rc.func_points = sorted(set(WenoConfig.func_points))
            rc.stencil_points, rc.function_stencil_dictionary = WenoConfig.generate_symbolic_function_points
            JS.generate_function_smoothness_indicators(rc)
            
            self.generate_alphas(rc, WenoConfig)
            self.generate_omegas(rc, WenoConfig)
            self.generate_reconstruction(rc, WenoConfig)
            self.reconstruction_classes[no] = rc

        # a = LeftReconstructionVariable('test')
        # a.function_stencil_dictionary = {}
        # t = DataSet('e')
        # fn = (DataSet('p')+DataSet('rhoE'))*DataSet('u0')
        # fn = DataSet('rhoE')*DataSet('u0')
        # for k in self.reconstruction_classes[0].func_points:
        #     loc = t.location
        #     loc[0] = loc[0] + k
        #     local_dict = {}
        #     for d in fn.atoms(DataSet):
        #         local_dict[d] = d.get_location_dataset(loc)
        #     a.function_stencil_dictionary[k] = fn.subs(local_dict)
        # pprint(a.function_stencil_dictionary)
        # self.reconstruction_classes[0].create_equations
        # a.create_function_equations(self.reconstruction_classes[0])
        # exit()
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
        fn.reconstructed_expression = stencil
        return

    def interpolate_reconstruction_variables(self, derivatives, kernel):
        for d in derivatives:
            pprint(d)
            for rv in d.reconstructions:
                if isinstance(rv, RightReconstructionVariable):
                    original_rv = self.reconstruction_classes[1]
                elif isinstance(rv, LeftReconstructionVariable):
                    original_rv = self.reconstruction_classes[0]
                else:
                    raise ValueError ("Reconstruction must be left or right")
                rv.update(original_rv)
                rv.add_evaluations_to_kernel(kernel)
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

    def convert_symbolic_to_dataset(self, symbolics, location, direction):

        symbols = symbolics.atoms(Symbol)
        dsets = [DataSet(s) for s in symbols]
        loc = dsets[0].location[:]
        loc[direction] = loc[direction] + location
        location_datasets = [d.get_location_dataset(loc) for d in dsets]
        substitutions = dict(zip(symbols, location_datasets))

        return symbolics.subs(substitutions)

class Characteristic(EigenSystem):
    def __init__(self,  eigenvalue, left_ev, right_ev):
        """ This holds the logic for the reconstruction procedure"""
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


    def pre_process(self, direction, derivatives, solution_vector, kernel):
        """
        """
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
               
        pre_process_equations += flatten(self.generate_equations_from_matrices(grid_LEV ,avg_LEV_values))
        
        # Transform the flux vector and the solution vector to characteristic space
        characteristic_flux_vector, CF_matrix = self.flux_vector_to_characteristic(derivatives, direction, averaged_suffix_name)
        characteristic_solution_vector, CS_matrix = self.solution_vector_to_characteristic(solution_vector, direction, averaged_suffix_name)
        pre_process_equations += flatten(self.generate_equations_from_matrices(CF_matrix, characteristic_flux_vector))
        pre_process_equations += flatten(self.generate_equations_from_matrices(CS_matrix, characteristic_solution_vector))     

        if hasattr(self, 'flux_split') and self.flux_split:
            print ("Characteristic flux splitting scheme")
            max_wavespeed_matrix, pre_process_equations = self.create_max_characteristic_wave_speed(pre_process_equations, direction)
            # self.split_fluxes(CF_matrix, CS_matrix, max_wavespeed_matrix)
            positive_flux = Rational(1,2)*(CF_matrix + max_wavespeed_matrix*CS_matrix)
            negative_flux = Rational(1,2)*(CF_matrix - max_wavespeed_matrix*CS_matrix)
            positive_flux.stencil_points = CF_matrix.stencil_points
            negative_flux.stencil_points = CF_matrix.stencil_points
            self.generate_right_reconstruction_variables(positive_flux, derivatives)
            self.generate_left_reconstruction_variables(negative_flux, derivatives)
        else:
            raise NotImplementedError("Only flux splitting is implemented in characteristic.")
        # self.pre_process_equations = pre_process_equations
        kernel.add_equation(pre_process_equations)
        return 

    def post_process(self, derivatives, kernel):
        post_process_equations = []
        averaged_suffix_name = self.averaged_suffix_name
        avg_REV_values = self.convert_matrix_to_grid_variable(self.right_eigen_vector[self.direction], averaged_suffix_name)
        grid_REV = self.generate_grid_variable_REV(self.direction, averaged_suffix_name) 
        post_process_equations += flatten(self.generate_equations_from_matrices(grid_REV ,avg_REV_values))
        reconstructed_characteristics = Matrix([d.evaluate_reconstruction for d in derivatives])
        reconstructed_flux = grid_REV*reconstructed_characteristics 
        reconstructed_work = [d.reconstruction_work for d in derivatives]
        post_process_equations += [Eq(x,y) for x,y in zip(reconstructed_work, reconstructed_flux)]
        kernel.add_equation(post_process_equations)
        
        return reconstructed_flux

    def generate_left_reconstruction_variables(self, flux, derivatives):
        if isinstance(flux, Matrix):
            stencil = flux.stencil_points
            for i in range(flux.shape[0]):
                rv = LeftReconstructionVariable('L_X%d_%d' % (self.direction, i))
                for p in sorted(set(self.reconstruction_classes[0].func_points)):
                    rv.function_stencil_dictionary[p] = flux[i,stencil.index(p)]
                derivatives[i].add_reconstruction_classes([rv])
        return

    def generate_right_reconstruction_variables(self, flux, derivatives):
        if isinstance(flux, Matrix):
            stencil = flux.stencil_points
            for i in range(flux.shape[0]):
                rv = RightReconstructionVariable('R_X%d_%d' % (self.direction, i))
                for p in sorted(set(self.reconstruction_classes[1].func_points)):
                    rv.function_stencil_dictionary[p] = flux[i,stencil.index(p)]
                derivatives[i].add_reconstruction_classes([rv])
        return 

    def solution_vector_to_characteristic(self, solution_vector, direction, name):
        stencil_points = sorted(list(set(self.reconstruction_classes[0].func_points + self.reconstruction_classes[1].func_points)))
        print "stencil points are: ", stencil_points
        solution_vector_stencil = zeros(len(solution_vector), len(stencil_points))
        CS_matrix = zeros(len(solution_vector), len(stencil_points))
        for j, val in enumerate(stencil_points): # j in fv stencil matrix
            for i, flux in enumerate(solution_vector):
                solution_vector_stencil[i,j] = self.increment_dataset(flux, direction, val)
                CS_matrix[i,j] = GridVariable('CS_%d%d' %(i,j))
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
        for j, val in enumerate(stencil_points): # j in fv stencil matrix
            for i, flux in enumerate(fv):
                flux_stencil[i,j] = self.increment_dataset(flux, direction, val)
                CF_matrix[i,j] = GridVariable('CF_%d%d' % (i,j))
        grid_LEV = self.generate_grid_variable_LEV(direction, name)
        characteristic_flux_stencil = grid_LEV*flux_stencil
        characteristic_flux_stencil.stencil_points = stencil_points
        CF_matrix.stencil_points = stencil_points
        return characteristic_flux_stencil, CF_matrix




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
        max_wave_speed = diag(*([max_wave_speed[i,i] for i in range(ev.shape[0])]))
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
        #print "Here "
        from .opensbliequations import *
        from .latex import *

        fname = './weno.tex'
        latex = LatexWriter()
        latex.open('./weno.tex')
        metadata = {"title": "Weno Flux equations", "author": "David", "institution": ""}
        latex.write_header(metadata)
    
        if isinstance(type_of_eq, SimulationEquations):
            eqs = flatten(type_of_eq.equations)
            grouped = self.group_by_direction(eqs)
            for key, derivatives in grouped.iteritems():
                pprint(key)
                for deriv in derivatives:
                    deriv.create_reconstruction_work_array(block)
                weno_kernel = Kernel(block)
                self.pre_process(key, derivatives, flatten(type_of_eq.time_advance_arrays), weno_kernel)
                pprint(len(weno_kernel.equations))
                self.interpolate_reconstruction_variables(derivatives, weno_kernel)
                
                reconstructed_flux = self.post_process(derivatives, weno_kernel)
                if key == 0:
                    weno_kernel.write_latex(latex)
                # self.post_process(derivatives)

            # pprint(type_of_eq.__dict__)
            latex.write_footer()
            latex.close()
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
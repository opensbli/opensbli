from sympy import *
from .opensblifunctions import *
from opensbli.core import *
from .grid import GridVariable
import time

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

class WenoSolutionType(object):
    def __init__(self, variable, left_right=None):
        self.variable = variable
        if left_right:
            self.reconstruction = left_right
        else:
            self.reconstruction = [False, True]
        self.reconstructedclass = {} 
        return
    def evaluate_interpolation(self, local_equations):
        symbols_list = []
        formulas_list = []
        reconstructed_symbols = {}
        # Left reconstruction
        if self.reconstruction[0]:
            left = self.reconstructedclass[-1]
            reconstructed_symbols[-1]= [left.symbolic_reconstructed]
            local_equations += left.all_equations
        # Right reconstruction
        if self.reconstruction[1]:
            right = self.reconstructedclass[1]
            reconstructed_symbols[1]= [right.symbolic_reconstructed]
            local_equations += right.all_equations
        return local_equations, reconstructed_symbols

class WenoConfig(object):
    def __init__(self, k, side):
        if side == -1:
            self.name, self.short_name, self.shift = 'left', 'L', 1
        elif side == 1:
            self.name, self.short_name, self.shift = 'right', 'R', 0
        self.k, self.side = k, side
        self.func_points = self.generate_left_right_points()
        # k passed explicitly as 2 sets of ENO coefficients are needed for the smoothness indicators
        self.c_rj = self.generate_eno_coefficients(k)
        self.c2_rj = self.generate_eno_coefficients(2*k-1)
        self.opt_weights = self.generate_optimal_weights()
        return

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

class Weno(Scheme):
    """ Main WENO class."""
    def __init__(self, order):
        """
        :returns: None
        """
        Scheme.__init__(self, "WenoDerivative", order)
        self.schemetype = "Spatial"
        self.k = int(0.5*(order+1))
        self.halotype = WenoHalos(order)
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        smooth_coeffs = self.generate_smoothness_coefficients()
        self.WenoConfigL = WenoConfig(self.k, -1)
        self.WenoConfigR = WenoConfig(self.k, 1)
        self.WenoConfigL.smooth_coeffs = smooth_coeffs
        self.WenoConfigR.smooth_coeffs = smooth_coeffs
        return

    def generate_weno(self, fn, direction, side, index):
        """ Generate the WENO scheme for a given function, direction and
        reconstruction point.

        :arg fn: The function to apply WENO to.
        :arg direction: Direction to apply WENO to.
        :arg side: Reconstruction point. Left: -1, Right: 1.
        :returns: The WENO reconstruction object.
        """
        self.direction, self.index = direction, index
        class reconstruction(object):
            pass
        # Store WENO attributes
        if side == -1:
            self.WC = self.WenoConfigL
        elif side == 1:
            self.WC = self.WenoConfigR
        weno = reconstruction()
        weno.all_equations = []
        # Generate the symbolic function points needed for the stencil
        weno.sym_func_points = self.generate_function_points(weno, fn)
        # Generate the weno coefficients
        weno.alpha = self.get_weno_coefficients(weno, fn)
        # Update the symbolic alphas used for reducing computations
        self.update_weno_equations(weno)
        # Create the final stencil
        weno.stencil = self.create_stencil(weno)
        weno.all_equations += [Eq(weno.symbolic_reconstructed, weno.stencil)]
        return weno

    def update_weno_equations(self, reconstruction):
        name, index = self.WC.short_name, self.index
        symbolic_list = []
        points_list = []
        for key, value in reconstruction.symbolc_points_dict.iteritems():
            symbolic_list += [value]
            points_list += [reconstruction.points_values[key]]
        reconstruction.all_equations += [Eq(a,b) for a,b in zip(symbolic_list, points_list)]
        # Beta equations
        reconstruction.all_equations += reconstruction.smoothness_equations
        # Alpha equations
        reconstruction.symbolic_alpha = [GridVariable('alpha_%s%d%d'%(name,index,i)) for i in range(0, self.k)]
        temp_zip = zip(reconstruction.symbolic_alpha, reconstruction.alpha)
        reconstruction.all_equations += [Eq(a,b) for a,b in temp_zip]
        reconstruction.sum_alpha =  GridVariable('Sigma_alpha_%s%d'%(name,index))
        reconstruction.symbolic_reconstructed = GridVariable('Recon_%s%d'%(name,index))
        reconstruction.all_equations += [Eq(reconstruction.sum_alpha, sum(reconstruction.symbolic_alpha))]
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

    def generate_function_points(self, reconstruction, expr):
        """ Indexes the function for a chosen direction and
        reconstruction point.

        arg: expr: The expression to apply WENO to.
        arg: direction: Direction to apply WENO.
        returns: all_fns: The function locations.
        """
        side, name, direction, index = self.WC.side, self.WC.name, self.direction, self.index
        reconstruction.symbolc_points_dict, reconstruction.points_values = {}, {}
        old_loc = [dset.location for dset in expr.atoms(DataSet)][0]
        for p in set(self.WC.func_points):
            if p>=0:
                reconstruction.symbolc_points_dict[p] = GridVariable('fn_%s%d_p%d'%(name,index,p))
            else:
                reconstruction.symbolc_points_dict[p] = GridVariable('fn_%s%d_m%d'%(name,index,abs(p)))
            # Get the current location of the function
            loc = old_loc[:]
            # Increment the location for each dataset based on the direction we are applying WENO to
            loc[direction] = old_loc[direction] + p
            for dset in expr.atoms(DataSet):
                expr = expr.replace(dset, dset.get_location_dataset(loc))
            reconstruction.points_values[p] = expr
        all_fns = [reconstruction.symbolc_points_dict[x] for x in self.WC.func_points]
        return all_fns

    def generate_smoothness_points(self, fn, shift, reconstruction):
        """ Creates the shifted function locations for the calculation of
        the smoothness indicator.

        arg: fn: The function to apply WENO to.
        arg: direction: Direction to apply WENO.
        arg: shift: The shift indices for the smoothness sum calculation.
        returns: fn_product: The product of the two function locations.
        """
        all_fns = []
        for p in shift:
            # if p in reconstruction.symbolc_points_dict.keys():
            all_fns += [reconstruction.symbolc_points_dict[p]]
        fn_product = all_fns[0]*all_fns[1]
        return fn_product

    def get_weno_coefficients(self, reconstruction, fn):
        """ Returns the WENO coefficients for a given order.

        arg: reconstruction: The WENO reconstruction object.
        arg: fn: The function to apply WENO to.
        arg: direction: Direction to apply WENO.
        returns: None
        """
        eps = Float(1e-16)
        smooth_coeffs, opt_weights, side = self.WC.smooth_coeffs, self.WC.opt_weights, self.WC.side
        k, name, shift, direction, index = self.k, self.WC.short_name, self.WC.shift, self.direction, self.index
        if k != 3:
            raise ValueError("WENO-Z is only implemented for k=3.") 
        smoothness_symbols = [GridVariable('beta_%s%d_%d'%(name, index, r)) for r in range(self.k)]
        smoothness_equations = []
        weno_z_symbol = GridVariable('tau_N_%s%d'%(name, index))
        alpha = []
        # Compute the smoothness indicator and alpha
        for r in range(k):
            smoothness_indicator = 0
            # Grid variable_format
            variable = smoothness_symbols[r]
            for m in range(0, 2*k-2): # WARNING: change to (0, 2*k -1) if this is now broken.
                for n in range(0, 2*k-2):
                    beta = smooth_coeffs.get((r,m,n))
                    if beta != None:
                        shift_indices = [-r+m+shift, -r+n+shift]
                        func_product = self.generate_smoothness_points(fn, shift_indices, reconstruction)
                        smoothness_indicator += beta*func_product
            smoothness_equations += [Eq(variable,smoothness_indicator)]
            alpha.append(opt_weights.get((0,r))*(Float(1.0) +weno_z_symbol/(eps +variable)))
        smoothness_equations += [Eq(weno_z_symbol, abs(smoothness_symbols[-1]- smoothness_symbols[0]))]
        reconstruction.smoothness_equations = smoothness_equations
        return alpha

    def create_stencil(self, reconstruction):
        """ Computes the linear sum of WENO & ENO coefficients & the function locations
         to create the final stencil.

         arg: reconstruction: The WENO reconstruction object.
         returns: stencil: The final WENO stencil.
         """
        stencil = 0
        k = self.k
        c_rj = self.WC.c_rj
        for r in range(k):
            eno_expansion = [c_rj.get((r,j))*reconstruction.sym_func_points[r*k+j] for j in range(k)]
            stencil += reconstruction.symbolic_alpha[r]*sum(eno_expansion)/reconstruction.sum_alpha
        return stencil

    def update_WenoSolutionType(self, required_interpolations):
        for no,interpolation in enumerate(required_interpolations):
            direction = interpolation.direction_index
            fn = interpolation.variable
            if interpolation.reconstruction[0]:
                side = -1
            if interpolation.reconstruction[1]:
                side = 1
            interpolation.sides = side
            interpolation.reconstructedclass[side] = self.generate_weno(fn, direction, side, no)
        return

class EigenSystem(object):
    def __init__(self, eigenvalue, left_ev, right_ev):
        if not isinstance(left_ev, dict) or not isinstance(right_ev, dict) or not isinstance(eigenvalue, dict):
            raise ValueError("Eigen values and eigen vectors should be in dictionary format")
        self.eigen_value = eigenvalue
        self.left_eigen_vector = left_ev
        self.right_eigen_vector = right_ev
        return

class Characteristic(EigenSystem):
    def __init__(self,  eigenvalue, left_ev, right_ev):
        """ This holds the logic for the reconstruction procedure"""
        self.is_vector_type = True
        self.leftRight = [True, True]
        EigenSystem.__init__(self, eigenvalue, left_ev, right_ev)
        return

    def generate_symbolic_LEV_REV(self):
        """
        Creates matrices containing the indexed LEV and REV placeholder symbols.
        """
        key = list(self.right_eigen_vector.keys())[0]
        matrixshape = self.right_eigen_vector[key].shape
        self.LEV_symbolic = zeros(*matrixshape)
        self.REV_symbolic = zeros(*matrixshape)
        for i in range(matrixshape[0]):
            for j in range(matrixshape[1]):
                self.LEV_symbolic[i,j] = GridVariable('LEV_%d%d'%(i,j))
                self.REV_symbolic[i,j] = GridVariable('REV_%d%d'%(i,j))
        return

    def convert_matrix_to_grid_variable(self, mat, name):
        """Converts the given matrix function to grid variable equivalent
        """
        syms = list(mat.atoms(Symbol))
        new_syms = [GridVariable('%s_%s'%(name,str(sym))) for sym in syms]
        substitutions = dict(zip(syms, new_syms))
        mat = mat.subs(substitutions)
        return mat

    def increment_dataset(self, eq, increment):
        """
        Increments a dataset by the given increment in the direction passed to the characteristic routine.
        """
        new_loc = self.base_location[:]
        new_loc[self.direction_index] += increment
        for dset in eq.atoms(DataSet):
            eq = eq.subs(dset, dset.get_location_dataset(new_loc))
        return eq

    def evaluate_eigen_values(self, formulas, required_symbols, side, position):
        """
        Evaluates rho, u0..u2, a for the eigenvalues for grid index i (left) and i+1 (right).
        """
        eqns = []
        if side == 'left':
            for eq in formulas:
                if eq.lhs in required_symbols:
                    left = eq.rhs
                    eqns += [Eq(GridVariable("%s_%s_%s"%(side, eq.lhs.args[0], position)), left)]
            eqns += [Eq(GridVariable("%s_%s_%s"%(side, self.rho.args[0], position)), self.rho)]
        elif side == 'right':
            for eq in formulas:
                if eq.lhs in required_symbols:
                    left = eq.rhs
                    right = self.increment_dataset(left, 1)
                    eqns += [Eq(GridVariable("%s_%s_%s"%(side, eq.lhs.args[0], position)), right)]
            eqns += [Eq(GridVariable("%s_%s_%s"%(side, self.rho.args[0], position)), self.increment_dataset(self.rho,1))]
        return eqns

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

    def simple_average_left_right(self, formulas, required_symbols, name):
        avg_equations = []
        # Loop over the formulas and check which ones are required for the eigenvalues
        for eq in formulas:
            if eq.lhs in required_symbols:
                # Increment by one to get left and right versions of each dataset (i and i+1)
                left = eq.rhs
                right = self.increment_dataset(left, 1)
                # Form the simple average 0.5*(f(i)+f(i+1))
                avg_equations += [Eq(GridVariable("%s_%s"%(name, eq.lhs.base)), Rational(1,2)*(left + right))]
        # Add the simple averaged rho
        left = self.rho
        right = self.increment_dataset(left, 1)
        avg_equations += [Eq(GridVariable("%s_%s"%(name, self.rho.base)), Rational(1,2)*(left + right))]
        return avg_equations

    def pre_process(self, direction, direction_index):
        """ 
        """
        self.direction = direction
        self.direction_index = direction_index
        pre_process_equations = []
        time_vector = (Matrix(self.vector_notation[CoordinateObject('t')]))
        # Save rho[0,0] for eigensystem evaluations
        self.rho = time_vector[0]
        # Location [0,0,0] in 3D
        self.base_location  = self.rho.location
        # Function versions of the symbols in the eigensystems
        required_ev_symbols = set([DataSet('a')] + [DataSet('rho')] + [DataSet('u%d' % i) for i in range(self.ndim)])
        # List of Euler formulas
        euler_formulas = self.required_formulas
        ## Eigenvalue name placeholders:
        #WARNING: Need to update the pre_process_equations with all of the correct placeholder EV names and equate to symbolic versions
        self.eigenvalue_names = {}
        names = ['leftm2', 'leftm1', 'LR', 'rightp1', 'rightp2', 'rightp3']
        for index, name in enumerate(names):
            self.eigenvalue_names[index-2] = diag(*[GridVariable('%s_lambda_%d' % (name, i)) for i in range(self.ndim+2)])
        # Evaluate left/right EV for m1, m2, p1, p2, p3
        evaluated_eigenvalue_quantities = []
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'left', 'm1')   
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'left', 'm2')
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'right', 'p1')
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'right', 'p2')
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'right', 'p3')
        # Perform simple average in rho, u, a
        name = 'LR'
        pre_processed_eqns = self.simple_average_left_right(euler_formulas, required_ev_symbols, name)
        # Convert eigenvalue matrix terms to LR_u0, LR_u0 + LR_a .. averaged symbols. 
        self.avg_eigen_values = self.convert_matrix_to_grid_variable(self.eigen_value[self.direction], name)
        # Create matrices of indexed LEV, REV placeholders
        self.generate_symbolic_LEV_REV()
        # Convert the LEV/REV value matrices to averaged LR symbols
        LEV = self.convert_matrix_to_grid_variable(self.left_eigen_vector[self.direction], name)
        REV = self.convert_matrix_to_grid_variable(self.right_eigen_vector[self.direction], name)
        # Create the equality between the symbols and the averaged values
        LEV_eqns = []
        REV_eqns = []
        for index in range(len(list(LEV))):
            LEV_eqns += [Eq(self.LEV_symbolic[index], LEV[index])]
            REV_eqns += [Eq(self.REV_symbolic[index], REV[index])]
            
        self.pre_process_equations = pre_processed_eqns
        self.pre_process_equations += LEV_eqns
        self.pre_process_equations += REV_eqns
        self.pre_process_equations += evaluated_eigenvalue_quantities
        # Interp multiplies the source vector by the LEV matrix elements, these are then passed to the update_weno solution routine
        # and then passed into the post processing part of the decomposition
        required_interpolations = self.interp_functions(direction)
        return required_interpolations

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
    def __init__(self,  eigenvalue, left_ev, right_ev, order, ndim):
        Characteristic.__init__(self, eigenvalue, left_ev, right_ev)
        Weno.__init__(self, order)
        self.ndim = ndim
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

class ScalarLocalLFScheme(Weno):
    def __init__(self, order, eqns, speeds, ndim):
        Weno.__init__(self, order)
        self.ndim = ndim
        self.speed = speeds
        self.grouped_eqns = self.group_by_direction(eqns)
        self.t = CoordinateObject('t')
        self.vector_notation = {}
        self.vector_notation[self.t] = self.get_time_derivative(eqns[0])
        self.get_space_derivatives()
        return

    def discretise(self, type_of_eq, block):
        ## Add the calling of functions for pre/post process and WENO in here
        return

    def get_time_derivative(self, eqns):
        """
        Get the time derivatives to add to the vector notation dictionary.
        """
        time_deriv = [deriv.args[0] for deriv in eqns.atoms(TemporalDerivative)]
        return time_deriv
    def get_space_derivatives(self):
        """
        Add space derivatives for each direction in the vector_notation dictionary.
        """
        coordinate_symbol = "x"
        cart = CoordinateObject('%s_i'%(coordinate_symbol))
        coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(self.ndim)]
        for index in coordinates:
            self.vector_notation[index] = []
        for i in range(self.ndim):
            if coordinates[i] in self.vector_notation:
                self.vector_notation[coordinates[i]] = [deriv.args[0] for deriv in self.grouped_eqns[i]]
        return 

    def group_by_direction(self, eqs):
        """
        Group the derivatives by direction.
        """
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

    def pre_process(self, direction, direction_index):
        """
        Find the Lax Friedrich fluxes.
        """
        # Direction x0, x1, x2 index 0, 1, 2
        self.direction_index = direction_index
        # Spatial and time vectors
        spatial_flux_vec = self.vector_notation[direction]
        time_vector = self.vector_notation[self.t]
        # Construct the fluxes
        self.fplus = Rational(1,2)*(Matrix(spatial_flux_vec) + Abs(self.speed[direction])*Matrix(time_vector))
        self.fminus = Rational(1,2)*(Matrix(spatial_flux_vec) - Abs(self.speed[direction])*Matrix(time_vector))
        required_interpolations = []
        # Required reconstruction
        leftRight = [False, True]
        # Perform the WENO procedure to each of the fluxes
        for flux in self.fplus:
            interpolation = WenoSolutionType(flux,leftRight)
            interpolation.direction = direction
            interpolation.direction_index = self.direction_index
            required_interpolations += [interpolation]
        leftRight = [True, False]
        for flux in self.fminus:
            interpolation = WenoSolutionType(flux,leftRight)
            interpolation.direction = direction
            interpolation.direction_index = self.direction_index
            required_interpolations += [interpolation]
        return required_interpolations

    def post_process(self,interpolated):
        self.post_process_equations = []
        temp_dictionary = {}
        for value in interpolated:
            self.post_process_equations, reconstructed_symbols = value.evaluate_interpolation(self.post_process_equations)
            temp_dictionary[value.variable] = reconstructed_symbols
        # The new naming uplus will be minus
        self.right_interpolated = []
        self.left_interpolated = []
        for val in self.fplus:
            self.right_interpolated += temp_dictionary[val][1]
        for val in self.fminus:
            self.left_interpolated += temp_dictionary[val][-1]
        final_flux  = (Matrix(self.right_interpolated) + Matrix(self.left_interpolated))
        print "final flux"
        pprint(final_flux)
        print "equations"
        pprint(self.post_process_equations)
        return self.post_process_equations,final_flux
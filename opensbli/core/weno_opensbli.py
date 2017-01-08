from sympy import *
from .opensblifunctions import *
from opensbli.core import *
import pyweno.symbolic as symbolic
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

class Weno(Scheme):
    """ The spatial derivatives of an arbitrary function 'F'
    on the numerical grid with the provided spatial scheme.

    For a wall boundary condition this will have a dependency on the grid range. """
    def __init__(self, order):
        """ Initialise the spatial derivative, which gives the equations
        of spatial Derivatives for the various combinations of the spatial scheme and order of accuracy.

        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :returns: None
        """
        Scheme.__init__(self, "WenoDerivative", order)
        self.schemetype = "Spatial"
        self.k = int(0.5*(order+1))
        self.halotype = WenoHalos(order)
        self.func_points = {}
        self.generate_left_right_points()
        pprint(self.func_points)
        return

    def generate_weno(self, fn, direction, side, index = None):
        """ Generate the WENO scheme for a given function, direction and
        reconstruction point.

        :arg fn: The function to apply WENO to.
        :arg direction: Direction to apply WENO to.
        :arg side: Reconstruction point. Left: -1, Right: 1.
        :returns: The WENO reconstruction object.
        Changing this to a list of equations
        """
        class reconstruction(object):
            pass
        if index:
            number = index
        else:
            number = 0
        weno = reconstruction()
        weno.all_equations = []
        pprint(fn)
        # Generate the function points needed for the stencil
        weno.func_points = self.generate_function_points(weno, fn, direction, side, number)
        # ENO coefficients
        weno.eno_coeffs = self.get_eno_coefficients(side)
        # Optimal coefficients for WENO
        weno.optimal_coeffs = self.get_optimal_coefficients(side)
        # Generate the weno coefficients
        weno.alpha = self.get_weno_coefficients(weno, fn, direction, side,number)
        # Update the symbolic alphas used for reducing computations
        # Create the final stencil
        self.update_weno_equations(weno, side, number)
        weno.stencil = self.create_stencil(weno)
        weno.all_equations += [Eq(weno.symbolic_reconstructed, weno.stencil)]
        return weno

    def update_weno_equations(self, reclass, side, number):
        if side == -1:
            name = 'L'
        elif side == 1:
            name = 'R' # This can be moved to the RECLASS LATER
        # First the equations for fn points
        ### EDIT THIS FUNCTION LATER
        symbolic_list = []
        points_list = []
        for key, value in reclass.symbolc_points_dict.iteritems():
            symbolic_list += [value]
            points_list += [reclass.points_values[key]]
        pprint(symbolic_list)
        pprint(points_list)
        reclass.all_equations += [Eq(a,b) for a,b in zip(symbolic_list, points_list)]
        # Beta equations
        reclass.all_equations += reclass.smoothness_equations
        # Alpha equations
        reclass.symbolic_alpha = [GridVariable('alpha_%s%d%d'%(name,number,i)) for i in range(0, self.k)]
        temp_zip = zip(reclass.symbolic_alpha, reclass.alpha)
        reclass.all_equations += [Eq(a,b) for a,b in temp_zip]
        reclass.sum_alpha =  GridVariable('Sigma_alpha_%s%d'%(name,number))
        reclass.symbolic_reconstructed = GridVariable('Recon_%s%d'%(name,number))
        reclass.all_equations += [Eq(reclass.sum_alpha, sum(reclass.symbolic_alpha))]
        return

    def generate_left_right_points(self):
        """ Populate the function evaluation points for left and right
        reconstruction.
        args: None
        returns: None
        """
        k = self.k
        r_values = [i for i in range(k)]
        fn_points = []
        # Left points
        side = -1
        for r in r_values:
            for j in range(k):
                fn_points.append(-r + j+1)
        self.func_points[-1] = fn_points
        fn_points = []
        # Right points
        side = 1
        for r in r_values:
            for j in range(k):
                fn_points.append((-r+j))
        self.func_points[1] = fn_points
        return

    def generate_function_points(self, reclass, fn, direction, side, number):
        """ Indexes the function for a chosen direction and
        reconstruction point.

        arg: fn: The function to apply WENO to.
        arg: direction: Direction to apply WENO.
        arg: side: Reconstruction point. Left: -1, Right: 1.
        returns: all_fns: The function locations.
        """
        all_indices = []
        wrt = direction
        # Save the functions as temporary grid variables
        if side == -1:
            name = 'L'
        elif side == 1:
            name = 'R'
        reclass.symbolc_points_dict = {}
        reclass.points_values = {}

        for p in set(self.func_points[side]):
            if p>=0:
                reclass.symbolc_points_dict[p] = GridVariable('fn_%s%d_p%d'%(name,number,p))
            else:
                reclass.symbolc_points_dict[p] = GridVariable('fn_%s%d_m%d'%(name,number,abs(p)))
            expr = fn
            # Get the current location of the function
            old_loc = [dset.location for dset in expr.atoms(DataSet)][0]
            loc = old_loc[:]
            # Increment the location based on the direction we are applying WENO to
            loc[direction] = old_loc[direction] + p
            for derivative in expr.args:
                base_func = derivative.args[-1]
                updated_derivative = derivative.replace(base_func, base_func.get_location_dataset(loc))
                expr = expr.replace(derivative, updated_derivative)
            reclass.points_values[p] = expr
        # pprint(reclass.symbolc_points_dict[p])
        # pprint(reclass.points_values)
        all_fns = [reclass.symbolc_points_dict[x] for x in self.func_points[side]]
        pprint(all_fns)
        return all_fns

    def generate_smoothness_points(self, fn, direction, shift, side, reclass):
        """ Creates the shifted function locations for the calculation of
        the smoothness indicator.

        arg: fn: The function to apply WENO to.
        arg: direction: Direction to apply WENO.
        arg: shift: The shift indices for the smoothness sum calculation.
        returns: fn_product: The product of the two function locations.
        """
        all_indices = []
        wrt = direction
        all_fns = []
         
        for p in shift:
            if p in reclass.symbolc_points_dict.keys():
                all_fns += [reclass.symbolc_points_dict[p]]
        fn_product = all_fns[0]*all_fns[1]
        return fn_product

    def get_eno_coefficients(self, side):
        """ Returns the ENO coefficients given the left/right
        reconstruction point.

        arg: side: Reconstruction point. Left: -1, Right: 1.
        returns: coeffs: Dictionary of ENO coefficients.
        """
        k = self.k
        coeffs = {}
        pyweno_output = symbolic.reconstruction_coefficients(k,[side])
        r_indices = [i for i in range(k)]
        # Store the ENO coefficients
        for r in r_indices:
            for j in range(k):
                index = (0,r,j)
                coeffs[index] = pyweno_output.get(index)
        return coeffs

    def get_optimal_coefficients(self, side):
        """ Returns the optimal WENO coefficients given the left/right
        reconstruction point.

        arg: side: Reconstruction point. Left: -1, Right: 1.
        returns: opt_coeffs: Dictionary of optimal WENO coefficients.
        """
        k = self.k
        opt_coeffs = {}
        pyweno_output = symbolic.optimal_weights(k,[side])
        # Store the optimal WENO coefficients
        for r in range(k):
            index = (0,r)
            opt_coeffs[index] = pyweno_output[0].get(index)
        return opt_coeffs

    def get_weno_coefficients(self, reconstruction, fn, direction, side, number):
        """ Returns the WENO coefficients for a given order.

        arg: reconstruction: The WENO reconstruction object.
        arg: fn: The function to apply WENO to.
        arg: direction: Direction to apply WENO.
        returns: None
        """
        # Constants for calculation of WENO coefficients
        eps = Float(1e-16)
        p = 2.0
        k = self.k
        # Smoothness coefficients
        smooth_coeffs = symbolic.jiang_shu_smoothness_coefficients(k)

        opt_coeffs = reconstruction.optimal_coeffs
        if side == -1:
            name = 'L'
        elif side == 1:
            name = 'R'

        smoothness_symbols = [GridVariable('beta_%s%d_%d'%(name, number, r)) for r in range(self.k)]

        smoothness_equations = []
        #### Add conditional here that WENO-Z only works for self.k == 3, otherwise use JS formulation without tau parameter
        weno_z_symbol = GridVariable('tau_N_%s%d'%(name, number))
        alpha = []
        # Compute the smoothness indicator and alpha

        for r in range(k):
            smoothness_indicator = 0
            # Grid variable_format
            variable = smoothness_symbols[r]
            for m in range(0, 2*k-1):
                for n in range(0, 2*k-1):
                    beta = smooth_coeffs.get((r,m,n))
                    if beta != None:
                        if side == 1:
                            shift = [-r+m, -r+n]
                        elif side == -1:
                            shift = [-r+m+1, -r+n+1]
                        func_product = self.generate_smoothness_points(fn, direction, shift, side, reconstruction)
                        smoothness_indicator += beta*func_product
            smoothness_equations += [Eq(variable,smoothness_indicator)]
            alpha.append(opt_coeffs.get((0,r))*(Float(1.0) +weno_z_symbol/(eps +variable)))

            self.alpha = alpha
        smoothness_equations += [Eq(weno_z_symbol, abs(smoothness_symbols[-1]- smoothness_symbols[0]))]
        reconstruction.smoothness_equations = smoothness_equations
        return  alpha

    def create_stencil(self, reconstruction):
        """ Computes the linear sum of WENO & ENO coefficients & the function locations
         to create the final stencil.

         arg: reconstruction: The WENO reconstruction object.
         returns: stencil: The final WENO stencil.
         """
        stencil = 0
        k = self.k
        for r in range(k):
            eno_expansion = []
            for j in range(k):
                eno_expansion.append(reconstruction.eno_coeffs.get((0,r,j))*reconstruction.func_points[r*k+j])
            stencil += reconstruction.symbolic_alpha[r]*sum(eno_expansion)/reconstruction.sum_alpha
        return stencil

    def update_WenoSolutionType(self, required_interpolations):
        for no,interpolation in enumerate(required_interpolations):
            direction = interpolation.direction_index
            fn = interpolation.variable
            sides = []
            if interpolation.reconstruction[0]:
                sides += [-1]
            if interpolation.reconstruction[1]:
                sides += [1]
            interpolation.sides = sides
            for side in sides:
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
        """ Derivatives should be of size of number of equations in flux vector format
        Characteristic decomposition evaluation requires
            # multiply with left Ev and interpolate them
        """
        self.direction = direction
        self.direction_index = direction_index
        pre_process_equations = []
        time_vector = (Matrix(self.vector_notation[CoordinateObject('t')]))
        # Save rho[0,0] for eigensystem evaluations
        self.rho = time_vector[0]
        # Location [0,0,0] in 3D
        self.base_location  = self.rho.location
        conservative_vars_base = [Symbol(str(l.base)) for l in time_vector]
        pprint(conservative_vars_base)
        # Solution vector at grid index i
        left = list(time_vector)
        # Solution vector at grid index i+1
        right = [self.increment_dataset(dset, 1) for dset in left]
        # Function versions of the symbols in the eigensystems
        required_ev_symbols = set([DataSet('a')] + [DataSet('rho')] + [DataSet('u%d' % i) for i in range(self.ndim)])
        # List of Euler formulas
        name = 'LR'
        euler_formulas = self.required_formulas
        # Perform simple average in rho, u, a
        pre_processed_eqns = self.simple_average_left_right(euler_formulas, required_ev_symbols, name)
        print "Simple averaged equations for eigenvalues are: "
        pprint(pre_processed_eqns)

        ######
        evaluated_eigenvalue_quantities = []
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'left', 'm1')   
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'left', 'm2')
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'right', 'p1')
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'right', 'p2')
        evaluated_eigenvalue_quantities += self.evaluate_eigen_values(euler_formulas, required_ev_symbols, 'right', 'p3')
        pprint(evaluated_eigenvalue_quantities)
        exit()


        # avg_eigen_name = 'LR'
        # self.eigenvalues_names = {}
        # # As the name is changed for to simple LR_*, convert the EV to this name
        # pre_process_equations += self.eigen_value_evaluation_eq(self.convert_to_grid_var_matrix(
        #     self.eigen_value[direction], name), avg_eigen_name)
        # pprint(pre_process_equations)
        # exit()
        ## Done up to here, need to create dictionary of diagonal matrices, one to store evs for each of the points in weno
        # at fifth order this is -2, -1, 0, 1, 2, 3, make it general for the order of WENO

        # Store the Eigen values as they are reused in the flux_evaluation
        # self.eigenvalues_names[0] = self.eigenvalues_symbolic

        # # Eigen values for i and i+1
        # # Substituting in left for m1, m2, right for p1, p2, p3 ()
        # pre_process_equations += self.evaluate_eigen_values(left_subs_dict, "left")
        # pre_process_equations += self.evaluate_eigen_values(left_subs_dict, "leftm2")
        # pre_process_equations += self.evaluate_eigen_values(right_subs_dict, "rightp1")
        # pre_process_equations += self.evaluate_eigen_values(right_subs_dict, "rightp2")
        # pre_process_equations += self.evaluate_eigen_values(right_subs_dict, "rightp3")
        # ### This is all naming stuff, leftm_2_lamda_0 = leftm2_symbol version etc
        # pre_process_equations += self.eigen_value_evaluation_eq(self.convert_to_grid_var_matrix(
        #     self.eigen_value[direction], "left"), "left")
        # self.eigenvalues_names[-1] = self.eigenvalues_symbolic
        # pre_process_equations += self.eigen_value_evaluation_eq(self.convert_to_grid_var_matrix(
        #     self.eigen_value[direction], "rightp1"), "rightp1")
        # self.eigenvalues_names[1] = self.eigenvalues_symbolic

        # pre_process_equations += self.eigen_value_evaluation_eq(self.convert_to_grid_var_matrix(
        #     self.eigen_value[direction], "rightp2"), "rightp2")
        # self.eigenvalues_names[2] = self.eigenvalues_symbolic

        # pre_process_equations += self.eigen_value_evaluation_eq(self.convert_to_grid_var_matrix(
        #     self.eigen_value[direction], "rightp3"), "rightp3")
        # self.eigenvalues_names[3] = self.eigenvalues_symbolic

        # pre_process_equations += self.eigen_value_evaluation_eq(self.convert_to_grid_var_matrix(
        #     self.eigen_value[direction], "leftm2"), "leftm2")
        # self.eigenvalues_names[-2] = self.eigenvalues_symbolic

        # This is calculating the LEV components and adding to pre_process
        # Left Eigen vector  matrix
        pre_process_equations += self.left_eigen_vectors_evaluation_eq(self.convert_to_grid_var_matrix(
            self.left_eigen_vector[direction], name))
        pprint(pre_process_equations)
        exit()
        # Calculating all of the REV components, in terms of LR_variable i.e. LEV31 = -sqrt2*(LR_a + 2*LR_u0)/(2*LR_a * LR_rho) <-------- this is where we
        # need an expression for rho, as it is in the LEV/REV matrices, are there any other terms I am missing in the new formulation of the eigensystems? 
        # LR_a for example is currently given by 'a' formula in terms of the simple averaged LR_a = 0.5*(formula for a at i, formula for a at i+1)
        # Right Eigen value matrix
        pre_process_equations += self.right_eigen_vectors_evaluation_eq(self.convert_to_grid_var_matrix(
            self.right_eigen_vector[direction], name))

        self.pre_process_equations = pre_process_equations

        ## Here the time vector (solution vector (rho, rhou0, rhou1, rhou2, rhoE) is passed to interp functions with the direction to interpolate by)
        # Need to see what interp_functions does. Self.direction also set to key, where is this used? required interpolations are what is returned by 
        # the pre_proecss function, all of the pre_process equations are stored to self. 

        # Interp multiplies the source vector by the LEV matrix elements, these are then passed to the update_weno solution routine
        # and then passed into the post processing part of the decomposition

        #characteristic = self.left_eigen_vector_symbolic*time
        exit()
        required_interpolations = self.interp_functions(time, key)
        pprint(required_interpolations)
        exit()
        self.direction = key


        return

class GLFCharacteristic(Characteristic, Weno):
    """ This class contains the Global Lax-Fedrich scheme
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

        This is main calling function from opensbli equations.spatial_discretisation, which is called from block.discretise
        Do here the following
        a. Find all the Function of type(self.name)
        b. Add all the DataSetBases required to the type_of_eq
        c. Create Equations for the evaluation ro create Kernels of each function depending on the grid/block
            control parameters
        d. Set the range of evaluation of the DataSetBases
        e. Update the Descritised equations in type_of_eq by substituting the equations with respective
            work array or discretised formula

        """  
                # self.set_halos(block)
        # from .opensbliequations import *
        # if isinstance(type_of_eq, SimulationEquations):
        #     if block.sbli_rhs_discretisation:
        #         self.sbli_rhs_discretisation(type_of_eq, block)
        #         return self.required_constituent_relations
        # else:
        #     local_kernels, discretised_eq = self.genral_discretisation(type_of_eq.equations, block)
        #     if discretised_eq:
        #         raise ValueError("")
        #     else:
        #         pass
        #     return self.required_constituent_relations    
        return
    def get_time_derivative(self, eqns):
        time_deriv = []
        for deriv in eqns.atoms(TemporalDerivative):
            time_deriv.append(deriv.args[0])
        return time_deriv

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

    def interp_functions(self, time_vector, key):
        local_equations = []
        f = lambda x:x.subs(x, abs(x))
        center = self.eigenvalues_names[0].applyfunc(f)
        # Center is diag matrix contianing abs value symbols LR_lambda_0..3 or 4
        diagonal = []
        for i in range(center.shape[0]):
            for j in range(center.shape[1]):
                if i == j:
                    diagonal += [abs(center[i,j])]
        # Diagonal is a list containing the diagonal terms of center
        max_lambda_names = [GridVariable('max_lambda')]
        local_equations = [Eq(max_lambda_names[0], Max(*(diagonal)))]
        lambda_max_matrix = diag(*([max_lambda_names[0]]*center.shape[0]))
        flux_vector = Matrix(self.vector_notation[key])
        censerv_vector = Matrix(self.vector_notation[EinsteinTerm('t')])
        ### lambda_max_matrix is a diagonal matrix saying "max_lambda" 4 times
        ### flux vector is flux vector from euler equations, self.vector(spatial direction we are using "here it is key")
        ### censerv vector is the solution vector rho, rhou0, rhou1, rhou2, rhoE

        # Maximum_lambda hard code them
        ### Here multiply the LEV symbolic matrix (LEV00 etc) by the flux and solution vectors
        ### Transformation to characteristic space is done here THIS IS PART (c) of the WENO procedure 2.10
        max_lam_vec = []
        flux_conserve = self.left_eigen_vector_symbolic*flux_vector
        censerv_vector = self.left_eigen_vector_symbolic*Matrix(self.vector_notation[EinsteinTerm('t')])
        ## Positive and negative are used to store ? 
        positive = []
        negative = []
        ## Loop over the length of flux, i= 4-5
        ## l_EV, r_EV etc store the names left_lambda_0 etc given for the 5 points, this needs to be generalised above
        ## In the pre_proecss function 
        ## This is the part where alpha is multiplied for the flux splitting f+ and f- are calculated here into positive, negative arrays
        ## The flux here is constructed once for each row of the flux/soln vectors f = 0.5*(f(u)+alpha*u), 4 times
        for i in range(len(flux_vector)):
            l_EV = self.eigenvalues_names[-1][i,i]
            r_EV = self.eigenvalues_names[1][i,i]
            c_EV = self.eigenvalues_names[0][i,i]
            r2EV = self.eigenvalues_names[2][i,i]
            r3EV = self.eigenvalues_names[3][i,i]
            lm2EV = self.eigenvalues_names[-2][i,i]
            # The solution for if local lambda == 0
            max_lambda = Max(Abs(c_EV),Abs(r_EV), Abs(l_EV), Abs(r2EV),Abs(lm2EV))#,Abs(r3EV))
            pprint(max_lambda)
            exit()
            # max(abs(lr_lambda_0)), abs(left_lambda_0) ... for each point, in a list? 
            # Adds "max_lambda_0, max_lambda_1  ... to a list"
            max_lam_vec += [GridVariable('max_lambda_%d'%i)]
            # Equates the above grid variable max_lambda_0 to the max_lambda equation expression for this component of flux vector
            # Adds to the massive pre_process list 
            self.pre_process_equations += [Eq(max_lam_vec[-1],max_lambda )]
            # This part does f+ = 0.5*(u(i) + MAX(alpha?)*f(u(i)))
            # and f- = 0.5*(u(i) - MAX(alpha?)*f(u(i)))
            # They have already been converted into characteristic space by multiplying by LEV in this function
            # This process is done 4 times in 2D, once for each component of the soln/flux vectors.  
            positive += [Rational(1,2)*(flux_conserve[i] + max_lam_vec[-1]*censerv_vector[i])]
            negative += [Rational(1,2)*(flux_conserve[i] - max_lam_vec[-1]*censerv_vector[i])]
        self.max_lam_vec = max_lam_vec

        # The fluxes in characteristic form are stored to self.u_plus/minus
        self.u_plus = positive
        self.u_minus = negative
        # Interpolations required  are right for uplus and left for uminus
        ## At this point the FULL flux equations in both u_plus/minus are passed to WENO to do the interpolations, with their key
        # left is used for u_minus (false, true)
        # right used for u_plus (true, false)

        ## ACTUAL weno procedure is applied in update_WenoSolutionType, which is called on the pre_processed equations in the main calling computations part
        ## This is done after WenoSolutionType is called updating which flux terms need WENO and the side to apply WENO to
        ## All of the below is setting up what weno procedures are needed, the actual weno is applied by update_WenoSolutionType.
        required_interpolations = []
        leftRight = [False, True]
        for val in self.u_plus:
            temp = WenoSolutionType(val,leftRight)
            temp.direction = key
            required_interpolations += [temp]
        leftRight = [True, False]
        for val in self.u_minus:
            temp = WenoSolutionType(val,leftRight)
            temp.direction = key
            required_interpolations += [temp]

        ### The interpolations are then returned back to weno_local control


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
        # self.set_halos(block)
        # from .opensbliequations import *
        # if isinstance(type_of_eq, SimulationEquations):
        #     if block.sbli_rhs_discretisation:
        #         self.sbli_rhs_discretisation(type_of_eq, block)
        #         return self.required_constituent_relations
        # else:
        #     local_kernels, discretised_eq = self.genral_discretisation(type_of_eq.equations, block)
        #     if discretised_eq:
        #         raise ValueError("")
        #     else:
        #         pass
        #     return self.required_constituent_relations
        # return
        return

    def get_time_derivative(self, eqns):
        """
        Get the time derivatives to add to the vector notation dictionary.
        """
        time_deriv = []
        for deriv in eqns.atoms(TemporalDerivative):
            time_deriv.append(deriv.args[0])
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
        print "finished post_process"
        return self.post_process_equations,final_flux

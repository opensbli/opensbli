from sympy import *
from .opensblifunctions import CentralDerivative
from .opensblifunctions import WenoDerivative
from .opensblifunctions import TemporalDerivative
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
            print "in left_right"
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
            print "in right"
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
        pprint(weno.func_points)
        exit()
        # ENO coefficients
        weno.eno_coeffs = self.get_eno_coefficients(side)
        # Optimal coefficients for WENO
        weno.optimal_coeffs = self.get_optimal_coefficients(side)
        # Generate the weno coefficients
        weno.alpha = self.get_weno_coefficients(weno, fn, direction, side,number)
        pprint(weno.alpha)
        # Update the symbolic alphas used for reducing computations
        # Create the final stencil
        self.update_weno_equations(weno, side, number)
        weno.stencil = self.create_stencil(weno)
        pprint(weno.stencil)
        # exit()
        weno.all_equations += [Eq(weno.symbolic_reconstructed, weno.stencil)]
        pprint(weno.all_equations)
        return weno

    def update_weno_equations(self, reclass, side, number):
        if side == -1:
            name = 'L'
        elif side == 1:
            name = 'R' # This can be moved to the RECLASS LATER
        # First the equations for fn points
        symbolic_list = []
        points_list = []
        for key, value in reclass.symbolc_points_dict.iteritems():
            symbolic_list += [value]
            points_list += [reclass.points_values[key]]
        print "here"
        pprint(symbolic_list)
        pprint(points_list)
        reclass.all_equations += [Eq(a,b) for a,b in zip(symbolic_list, points_list)]
        # Beta equations
        reclass.all_equations += reclass.smoothness_equations
        pprint(reclass.smoothness_equations)
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
            pprint(p)
            if p>=0:
                reclass.symbolc_points_dict[p] = GridVariable('fn_%s%d_p%d'%(name,number,p))
                pprint(reclass.symbolc_points_dict[p])
            else:
                reclass.symbolc_points_dict[p] = GridVariable('fn_%s%d_m%d'%(name,number,abs(p)))
            expr = fn
            # Get the current location of the function
            old_loc = [dset.location for dset in expr.atoms(DataSet)][0]
            loc = old_loc[:]
            # Increment the location based on the direction we are applying WENO to
            loc[direction] = old_loc[direction] + 999 # 999 should be p when not testing (p starts at 0)

            #WARNING : I don't know how to update the locations for both terms in the flux fn, while retaining
            # the original function of 1/2 phi[0,0] + |c_0| phi[0,0]/2
            # e.g. (1/2 phi[0,0] + |c_0| phi[0,0]/2) -----> 1/2 phi[999,0] + |c_0| phi[999,0]/2
            # old location [0,0] is stored in old_loc, new location [999,0] stored in loc, expr needs to be updated
            # run wave.py to get output up until here
            pprint(srepr(expr))
            pprint(fn)
            pprint(expr.args)
            pprint(old_loc)
            pprint(loc)
            for func in expr.args:
                func.args[-1].replace(old_loc, loc)
                pprint(func.args)
            pprint(func)
            exit()
            expr = expr.subs(old_loc, loc)
            pprint(expr)
            exit()
            reclass.points_values[p] = expr

        pprint(reclass.points_values)
        all_fns = [reclass.symbolc_points_dict[x] for x in self.func_points[side]]
        pprint(all_fns)
        # exit()
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


class GLFCharacteristic(Characteristic, Weno):
    """ This class contains the Global Lax-Fedrich scheme
    """
    def __init__(self,  eigenvalue, left_ev, right_ev, order):
        Characteristic.__init__(self, eigenvalue, left_ev, right_ev)
        Weno.__init__(self, order)
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

class ScalarLocalLFScheme(Weno):

    def __init__(self, order, eqns, speeds, ndim):
        Weno.__init__(self, order)
        self.ndim = ndim
        self.speed = speeds
        self.grouped_eqns = self.group_by_direction(eqns)

        self.t = ConstantObject('t')
        self.vector_notation = {}
        self.vector_notation[self.t] = self.get_time_derivative(eqns[0])
        self.get_space_derivatives()
        pprint(self.vector_notation)
        return

    def discretise(self, type_of_eq, block):
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
        pprint(self.fplus)
        pprint(self.fminus)
        required_interpolations = []
        # Required reconstruction
        leftRight = [False, True]
        for flux in self.fplus:
            pprint(flux)
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
        pprint(temp_dictionary)
        exit()
        # The new naming uplus will be minus
        self.right_interpolated = []
        self.left_interpolated = []
        for val in self.fplus:
            self.right_interpolated += temp_dictionary[val][1]
        for val in self.fminus:
            self.left_interpolated += temp_dictionary[val][-1]
        final_flux  = (Matrix(self.right_interpolated) + Matrix(self.left_interpolated))
        return self.post_process_equations,final_flux
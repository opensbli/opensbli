
# from .scheme import Scheme

class WenoHalos(object):
    def __init__(self, order):
        # Check for the boundary types in the blocks and set the halo points
        #self.halos = [[-scheme.order, scheme.order] for dim in range(block.ndim)]
        return
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
        # self.derivative_direction = grid.indices
        # self.index_mapping = grid.mapped_indices
        Scheme.__init__(self, "WenoDerivative", order)
        self.schemetype = "Spatial"
        self.k = int(0.5*(order+1))
        self.halotype = WenoHalos(order)
        print "Inside Weno"
        return

class EigenSystem(object):
    def __init__(self, eigenvalue, left_ev, right_ev):
        if not isinstance(left_ev, dict) or not isinstance(right_ev, dict) or not isinstance(eigenvalue, dict):
            raise ValueError("Eigen values and eigen vectors should be in dictionary format")
        # TODO the above line should be made more explicit one condition for each input
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
        #TODO fix the halo points here
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

class ScalarLocalLFScheme(Weno):

    def __init__(self, order):
        self.is_vector_type = True
        Weno.__init__(self, order)
        return

    # def halos_required(self, k, ndim):
    #     halos = [(-k, k+1)for dim in range(ndim)]
    #     return halos
    # def set_fluxes_direction(self, eq_vector_form):
    #     self.vector_notation = eq_vector_form
    #     return
    # def set_scalar_speed(self, speed):
    #     self.speed = speed
    #     return
    # def pre_process(self, key):
    #     """ Find the lax fedrich fluxes"""
    #     spatial_flux_vec = self.vector_notation[key]
    #     time_vector = self.vector_notation[EinsteinTerm('t')]
    #     pprint(time_vector)
    #     pprint(Abs(self.speed[key]))
    #     # find the fluxes
    #     self.fplus = Rational(1,2)*(Matrix(spatial_flux_vec) + Abs(self.speed[key])*Matrix(time_vector))
    #     self.fminus = Rational(1,2)*(Matrix(spatial_flux_vec) - Abs(self.speed[key])*Matrix(time_vector))
    #     required_interpolations = []
    #     leftRight = [False, True]
    #     for val in self.fplus:
    #         temp = WenoSolutionType(val,leftRight)
    #         temp.direction = key
    #         required_interpolations += [temp]
    #     leftRight = [True, False]
    #     for val in self.fminus:
    #         temp = WenoSolutionType(val,leftRight)
    #         temp.direction = key
    #         required_interpolations += [temp]
    #     return required_interpolations
    # def post_process(self,interpolated):
    #     self.post_process_equations = []
    #     temp_dictionary = {}
    #     for val in interpolated:
    #         self.post_process_equations, reconstructed_symbols = val.evaluate_interpolation(self.post_process_equations)
    #         temp_dictionary[val.variable] = reconstructed_symbols
    #     #pprint(temp_dictionary)
    #     # The new naming uplus will be minus
    #     self.right_interpolated = []
    #     self.left_interpolated = []
    #     for val in self.fplus:
    #         self.right_interpolated += temp_dictionary[val][1]
    #     for val in self.fminus:
    #         self.left_interpolated += temp_dictionary[val][-1]
    #     final_flux  = (Matrix(self.right_interpolated) + Matrix(self.left_interpolated))
    #     return self.post_process_equations,final_flux
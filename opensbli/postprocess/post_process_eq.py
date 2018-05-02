from opensbli.core.parsing import EinsteinEquation
from sympy import Eq, zeros, flatten, Matrix, pprint, Function, S, Equality, Wild, WildFunction, srepr
from opensbli.initialisation.common import AfterSimulationEnds
from opensbli.core.opensbliequations import NonSimulationEquations, Discretisation, Solution, OpenSBLIEquation, DataSet
from opensbli.core.opensblifunctions import CentralDerivative, WenoDerivative
from sympy.tensor.array import MutableDenseNDimArray
from opensbli.core.opensbliobjects import CoordinateObject, DataObject
from opensbli.core.kernel import Kernel
from opensbli.core.latex import LatexWriter
from opensbli.core.bcs import BoundaryConditionBase, ModifyCentralDerivative, Carpenter


class DummyCarpenter(ModifyCentralDerivative, BoundaryConditionBase):
    """ Navier-Stokes specific boundary condition. Applies a no-slip viscous wall condition, 
    velocity components are zero on the wall. Temperature is fixed with a prescribed wall temperature, 
    given as the rhoE equation passed to this BC. A wall normal velocity is set based on the equation passed by the user.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg Eq equations: Equation for conservative variable rhoE to set on the wall with constant temperature.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""
    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'Dummy'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        return []

class PostProcessEquations(NonSimulationEquations, Discretisation, Solution):
    """Class for the simulation equations. This performs the discretisation of the equations.
    
    :param int order: priority in the algorithm if multiple simulation equations exitst
    """
 
    def __new__(cls, order=None, place=None, **kwargs):
        ret = super(PostProcessEquations, cls).__new__(cls)
        #ret = object.__new__(cls)
        if order:
            ret.order = order
        else:
            ret.order = 0
        if place:
            ret.algorithm_place = [place]
        else:
            # Default to after the simulation
            ret.algorithm_place = [AfterSimulationEnds()]
        ret.equations = []
        ret.kwargs = kwargs
        return ret

    def add_constituent_relations(cls, constituent_relations):
        """Convert the constituent relations passed to a dictionary for easier access in discretisation
        
        :param ConstituentRelations constituent_relations: the constituent relations class on the block"""

        cls.constituent_relations_dictionary = cls.convert_to_dictionary(constituent_relations.equations)
        return

    def create_residual_arrays(cls):
        """Creates the residual datasets for each of the simulation equations.
        
        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        for eq in flatten(cls.equations):
            eq.residual = eq.lhs
        return
    
    @property
    def solution_arrays(cls):
        return [eq.residual for eq in flatten(cls.equations)]

    @property
    def get_required_constituents(self):
        arrays = []
        for eq in flatten(self.equations):
            arrays += list(eq.atoms(DataSet))
        arrays = set(arrays)
        return arrays

    def spatial_discretisation(cls, block):
        """Apllies the spatial discretisation of the equations by calling the discretisation of each spatial 
        scheme provided on the block
        
        :param SimulationBlock block: the block on which the equations are solved
        :return: None """

        # Instantiate the solution class
        #(Solution, cls).__init__(cls)
        cls.solution = Solution()
        # Create the residual array for the equations
        cls.create_residual_arrays()

        spatialschemes = []
        # Get the schemes on the block
        schemes = block.discretisation_schemes
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation
        cls.constituent_evaluations = {}
        crs = block.get_constituent_equation_class
        cr_dictionary = {}
        for cr in crs:
            cr_dictionary.update(cr.get_relations_dictionary)
        cls.requires = {}
        for no, sc in enumerate(spatialschemes):
            cls.constituent_evaluations[sc] = schemes[sc].discretise(cls, block)
            for key, value in cls.constituent_evaluations[sc].iteritems():
                if key in cr_dictionary.keys():
                    if key in cls.constituent_relations_kernels:
                        cls.constituent_relations_kernels[key].merge_halo_range(value.halo_ranges)
                    else:
                        cls.constituent_relations_kernels[key] = value
                        cls.constituent_relations_kernels[key].add_equation(cr_dictionary[key])
                else:
                    # raise ValueError("Constituent relation is not found for %s"%key)
                    cls.requires[key] = value
        missing_CR_datasets = cls.get_required_constituents.difference(cls.constituent_relations_kernels.keys())
        for dset in missing_CR_datasets:
            # Evaluation of missing dataset is required
            if dset in cr_dictionary.keys():
                for kernel in cr_dictionary[dset].kernels:
                    cls.constituent_relations_kernels[kernel.equations[0].lhs] = kernel
        cls.Kernels = cls.sort_constituents + cls.Kernels
        cls.process_kernels(block)
        return

    def process_kernels(cls, block):
        """A function to update some dependant parameters of each kernel
        
        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        for key, kernel in cls.constituent_relations_kernels.iteritems():
            kernel.update_block_datasets(block)
        for kernel in cls.Kernels:
            kernel.update_block_datasets(block)
        return

    @property
    def sort_constituents(cls):
        """Sort the constituent relation kernels
        """
        input_order = []
        for a in cls.requires.keys():
            if isinstance(a, DataSet):
                input_order += [a.base]
            else:
                input_order += [a]

        dictionary = {}
        order_of_evaluation = []
        for key, value in cls.constituent_relations_kernels.iteritems():
            dictionary[key.base] = value
            order_of_evaluation += [key.base] 
        #order_of_evaluation = cls.sort_dictionary(input_order, dictionary)
        ordered_kernels = []
        for o in order_of_evaluation:
            ordered_kernels += [dictionary[o]]
        return ordered_kernels

    def sort_dictionary(cls, order, new_dictionary):
        """ Sort the evaluations based on the requirements of each term. For example, if we have
        the primitive variables p, u0, u1, and T, then the pressure p may depend on the velocity u0 and u1, and T may depend on p,
        so we need this be evaluate in the following order: u0, u1, p, T.


        :arg list order: The list of already sorted terms or Known terms
        :arg evaluations: The evaluation information, containing dependency information.
        :arg typef: The type of term to sort.
        :returns: A list of ordered terms.
        :rtype: list
        """
        dictionary = new_dictionary
        # reverse_dictionary = {}
        order = flatten(order + [a.base for a in flatten(cls.solution_arrays)])
        order = list(set(order))
        # store the length of order
        input_order = len(order)
        key_list = [key for key in dictionary.keys() if key not in order]
        requires_list = ([dictionary[key].required_data_sets for key in key_list])

        zipped = zip(key_list, requires_list)
        # Breaks after 1000 iterations
        iter_count = 0
        while key_list:
            iter_count = iter_count+1
            order += [x for (x, y) in zipped if all(req in order for req in y)]
            key_list = [key for key in dictionary.keys() if key not in order]
            requires_list = [dictionary[key].required_data_sets for key in key_list]
            zipped = zip(key_list, requires_list)
            if iter_count > 1000:
                print("Exiting because i cannot classify the following")
                print("Already sorted are")
                pprint(order)
                pprint([srepr(o) for o in order])
                # print("Trying to sort the required for")
                # pprint(evaluations[key].lhs)
                print("It requires")
                pprint([req for req in requires_list[0]])
                print("Sorted")
                pprint([(req, req in order) for req in requires_list[0]])
                raise ValueError("Exiting sort evaluations ")
        order = order[input_order:]
        return order

    def apply_boundary_conditions(cls, block):
        """ Used to fit in the abstraction, boundary conditions are applied for each equation class passed to block and discretised.
        Later once metrics are not part of the eqaution classes, this can be removed. """
        return

    def apply_interface_bc(cls, block, multiblock_descriptor):
        return

    @property
    def all_spatial_kernels(cls):
        return cls.sort_constituents + cls.Kernels

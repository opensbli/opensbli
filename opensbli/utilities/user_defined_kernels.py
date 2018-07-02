from opensbli.equation_types.opensbliequations import NonSimulationEquations, Discretisation, Solution
from opensbli.core.kernel import Kernel


class UserDefinedEquations(NonSimulationEquations, Discretisation, Solution):
    """User defined equations, this will not discretise the equations. No checking is performed.
    Just forms a kernel on the range and places the kernel in the algorithm place passed as an
    input to the class
    """

    def __new__(cls, **kwargs):
        ret = super(UserDefinedEquations, cls).__new__(cls)
        ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        ret._place = []
        return ret

    @property
    def algorithm_place(cls):
        return cls._place

    @algorithm_place.setter
    def algorithm_place(cls, place):
        cls._place += [place]
        return

    def spatial_discretisation(cls, block):
        """ Applies the spatial discretisation of the equations by calling the discretisation of each spatial scheme provided on the block

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """

        # Instantiate the solution class
        cls.solution = Solution()
        user_defined_kernel = Kernel(block)
        name = " ".join([p.__class__.__name__ for p in cls.algorithm_place])
        user_defined_kernel.set_computation_name("user kernel %s" % name)
        user_defined_kernel.set_grid_range(block)
        user_defined_kernel.add_equation(cls.equations)
        user_defined_kernel.update_block_datasets(block)
        cls.Kernels += [user_defined_kernel]
        return

    def apply_boundary_conditions(cls, block):
        return

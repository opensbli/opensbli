from opensbli.equation_types.opensbliequations import NonSimulationEquations, Discretisation, Solution
from opensbli.core.kernel import Kernel
from sympy import flatten, pprint, Equality
from opensbli.core.opensbliobjects import GroupedPiecewise


class UserDefinedEquations(NonSimulationEquations, Discretisation, Solution):
    """User defined equations. No checking is performed.
    Just forms a kernel on the range and places the kernel in the algorithm place passed as an
    input to the class. """

    def __new__(cls, **kwargs):
        ret = super(UserDefinedEquations, cls).__new__(cls)
        ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        ret._place = []
        ret.computation_name = None
        # Optional halo type
        ret.halos = None
        return ret

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        for eq in flatten(cls.equations):
            if isinstance(eq, Equality):
                evaluated = evaluated.union(eq.lhs_datasetbases)
            elif isinstance(eq, GroupedPiecewise):
                evaluated = evaluated.union(eq.lhs_datasetbases)
        return evaluated

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
        # cls.solution = Solution()

        # Discretize any derivatives in the equations
        spatialschemes = []
        # Get the schemes on the block
        schemes = block.discretisation_schemes
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation if any in constituent relations evaluation
        equations = cls.equations

        UDF_derivative_kernels = []
        evaluations = []

        for eq in flatten(equations):
            cls.equations = [eq]
            for sc in spatialschemes:
                evaluations.append(schemes[sc].discretise(cls, block))
            UDF_derivative_kernels.append(cls.Kernels[:])

            cls.Kernels = []
        cls.equations = equations

        # No discretisation required
        if not flatten(UDF_derivative_kernels):
            user_defined_kernel = Kernel(block)
            user_defined_kernel.set_computation_name("User kernel: %s" % str(cls.computation_name))
            user_defined_kernel.set_grid_range(block)
            # Evaluate the kernel into the halos as required
            if cls.halos:
                for direction in range(block.ndim):
                    user_defined_kernel.set_halo_range(direction, 0, cls.halos[direction][0])
                    user_defined_kernel.set_halo_range(direction, 1, cls.halos[direction][1])
            user_defined_kernel.add_equation(cls.equations)
            # user_defined_kernel.update_block_datasets(block)
            cls.Kernels += [user_defined_kernel]
        else:
            cls.Kernels += flatten(UDF_derivative_kernels)

        # Process the kernels to update parameters on the block
        cls.process_kernels(block)
        return

    def process_kernels(cls, block):
        """A function to update some dependant parameters of each kernel

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        for kernel in cls.Kernels:
            kernel.update_block_datasets(block)
        return

    def apply_boundary_conditions(cls, block):
        return

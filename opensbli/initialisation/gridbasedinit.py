from opensbli.core.opensbliequations import NonSimulationEquations, Discretisation, OpenSBLIEquation
from opensbli.core.kernel import Kernel
from .common import BeforeSimulationStarts


class GridBasedInitialisation(NonSimulationEquations):
    def __new__(cls, order=None, **kwargs):
        ret = super(GridBasedInitialisation, cls).__new__(cls)
        if order:  # Local order if multiple instances of the class are declared on the block
            ret.order = order
        else:
            ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        """A control parameter is needed for where to put these equations in the algorithm
        """
        ret.algorithm_place = [BeforeSimulationStarts()]
        return ret

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return "GridBasedInitialisation"

    def add_equations(cls, equation):
        # equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq.lhs, eq.rhs)
                # eq.set_vector(no)
                cls.equations += [eq]
        else:
            equation = OpenSBLIEquation(equation.lhs, equation.rhs)
            cls.equations += [equation]
        return

    def spatial_discretisation(cls, block):
        kernel1 = Kernel(block, computation_name="Grid_based_initialisation%d" % cls.order)
        kernel1.set_grid_range(block)
        # Checking
        # for eq in block.list_of_equation_classes:
        # from opensbli.core.metrics import *
        # if isinstance(eq, MetricsEquation):
        # for k in eq.Kernels:
        # print k.computation_name, k.halo_ranges, k.ranges
        schemes = block.discretisation_schemes
        for d in range(block.ndim):
            for sc in schemes:
                if schemes[sc].schemetype == "Spatial":
                    kernel1.set_halo_range(d, 0, schemes[sc].halotype)
                    kernel1.set_halo_range(d, 1, schemes[sc].halotype)
        kernel1.add_equation(cls.equations)
        kernel1.update_block_datasets(block)
        cls.Kernels = [kernel1]
        # Checking
        # for eq in block.list_of_equation_classes:
        # from opensbli.core.metrics import *
        # if isinstance(eq, MetricsEquation):
        # for k in eq.Kernels:
        # print k.computation_name, k.halo_ranges, k.ranges
        return

    def apply_boundary_conditions(cls, block):
        """No boundary conditions in the Initialisation
        """
        return


# class GridGeneration(Discretisation, NonSimulationEquations):
#     def __new__(cls, )

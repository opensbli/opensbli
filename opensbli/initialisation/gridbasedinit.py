from opensbli.equation_types.opensbliequations import NonSimulationEquations, OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.opensbliobjects import GroupedPiecewise
from sympy import Equality
from opensbli.code_generation.algorithm.common import BeforeSimulationStarts


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

    def add_equations(cls, equation):  # TODO: Check if this is required, nonsimulation equations are derived from discretisation class
        # equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            for no, eq in enumerate(equation):
                if isinstance(eq, Equality):
                    eq = OpenSBLIEq(eq.lhs, eq.rhs)
                    cls.equations += [eq]
                elif isinstance(eq, GroupedPiecewise):
                    cls.equations += [eq]
        else:
            if isinstance(equation, Equality):
                eq = OpenSBLIEq(equation.lhs, equation.rhs)
                cls.equations += [eq]
            elif isinstance(equation, GroupedPiecewise):
                cls.equations += [equation]
            else:
                raise TypeError("Provide the equations in a list")
            # equation = OpenSBLIEquation(equation.lhs, equation.rhs)
            # cls.equations += [equation]
        return

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        for eq in cls.equations:
            if isinstance(eq, Equality):
                evaluated = evaluated.union(eq.lhs_datasetbases)
            elif isinstance(eq, GroupedPiecewise):
                evaluated = evaluated.union(eq.lhs_datasetbases)
        return evaluated

    def spatial_discretisation(cls, block):
        kernel1 = Kernel(block, computation_name="Grid_based_initialisation%d" % cls.order)
        kernel1.set_grid_range(block)
        schemes = block.discretisation_schemes
        for d in range(block.ndim):
            for sc in schemes:
                if schemes[sc].schemetype == "Spatial":
                    kernel1.set_halo_range(d, 0, schemes[sc].halotype)
                    kernel1.set_halo_range(d, 1, schemes[sc].halotype)
        kernel1.add_equation(cls.equations)
        kernel1.update_block_datasets(block)
        cls.Kernels = [kernel1]
        return

    def apply_boundary_conditions(cls, block):
        """No boundary conditions in the Initialisation
        """
        return

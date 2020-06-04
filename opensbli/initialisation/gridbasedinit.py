"""@brief
   @author Satya Pramod Jammy
   @contributors David Lusher
   @details
"""

from opensbli.equation_types.opensbliequations import NonSimulationEquations
from opensbli.core.kernel import Kernel
from opensbli.core.opensbliobjects import GroupedPiecewise
from sympy import Equality, flatten
from opensbli.code_generation.algorithm.common import BeforeSimulationStarts
from opensbli.schemes.spatial.scheme import CentralHalos_defdec


class GridBasedInitialisation(NonSimulationEquations):
    def __new__(cls, order=None, **kwargs):
        ret = super(GridBasedInitialisation, cls).__new__(cls)
        if order:  # Local order if multiple instances of the class are declared on the block more for future work
            ret.order = order
        else:
            ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        # A control parameter is needed for where to put these equations in the algorithm
        ret.algorithm_place = [BeforeSimulationStarts()]
        return ret

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return "GridBasedInitialisation"

    @property
    def evaluated_datasets(cls):
        evaluated = set()
        for eq in flatten(cls.equations):
            if isinstance(eq, Equality):
                evaluated = evaluated.union(eq.lhs_datasetbases)
            elif isinstance(eq, GroupedPiecewise):
                evaluated = evaluated.union(eq.lhs_datasetbases)
        return evaluated

    def spatial_discretisation(cls, block):
        kernel = Kernel(block, computation_name="Grid_based_initialisation%d" % cls.order)
        kernel.set_grid_range(block)
        schemes = block.discretisation_schemes
        for d in range(block.ndim):
            # Initialize all five halos
            kernel.set_halo_range(d, 0, CentralHalos_defdec())
            kernel.set_halo_range(d, 1, CentralHalos_defdec())
        kernel.add_equation(cls.equations)
        kernel.update_block_datasets(block)
        cls.Kernels = [kernel]
        return

    def apply_boundary_conditions(cls, block):
        """No boundary conditions in the Initialisation currently, any logic needed should be implemented here
        """
        return

    def apply_interface_bc(cls, block, multiblock_descriptor):
        """ No BC for initialisation."""
        return



#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (c) see License file

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

"""@brief
   @author Satya Pramod Jammy
   @contributors David Lusher
   @details
"""

from opensbli.equation_types.opensbliequations import NonSimulationEquations, OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.opensbliobjects import GroupedPiecewise
from sympy import Equality
from opensbli.code_generation.algorithm.common import BeforeSimulationStarts


class GridBasedInitialisation(NonSimulationEquations):
    def __new__(cls, order=None, **kwargs):
        ret = super(GridBasedInitialisation, cls).__new__(cls)
        if order:  # Local order if multiple instances of the class are declared on the block more for future work
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
        """No boundary conditions in the Initialisation currently, any logic needed should be implemented here
        """
        return

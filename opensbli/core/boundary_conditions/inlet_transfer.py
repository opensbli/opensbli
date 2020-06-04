from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from sympy import flatten


class InletTransferBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Simple inlet boundary condition to copy all solution variable values from the left halos
    to the boundary plane.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'InletTransfer'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        if side != 0:
            raise ValueError("Only implemented this BC for inlet side 0.")
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        cons_vars = flatten(arrays)
        equations = self.create_boundary_equations(cons_vars, cons_vars, [(0, -1)])
        kernel.add_equation(equations)
        kernel.update_block_datasets(block)
        return kernel

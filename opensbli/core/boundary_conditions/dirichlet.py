from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.opensbliobjects import ConstantIndexed
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter


class DirichletBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Applies a constant value Dirichlet boundary condition.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg list equations: OpenSBLI equations to enforce on the boundary.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, equations, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'Dirichlet'
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def halos(self):
        return True

    def apply(self, arrays, block):
        direction, side = self.direction, self.side
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        # Dirichlet set on the boundary
        kernel.add_equation(self.equations)
        # Change ranges if using split BC
        if isinstance(kernel.halo_ranges, ConstantIndexed):
            # Manually set Dirichlet into the halo range of this side
            halo_object = kernel.halo_ranges
            halo_object._value[2*direction + side] = halos[direction][side]
        else:  # Not using split BC, halos should be updated
            kernel.halo_ranges[direction][side] = block.boundary_halos[direction][side]
        kernel.update_block_datasets(block)
        return kernel

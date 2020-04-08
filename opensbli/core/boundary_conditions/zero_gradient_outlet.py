from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase
from sympy import flatten


class ZeroGradientOutletBC(BoundaryConditionBase):
    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'ZeroGradientOutlet'
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        n_halos = abs(halos[self.direction][self.side])
        from_side, to_side = self.set_side_factor()
        transfer_indices = [(from_side*i, to_side*i) for i in range(1, n_halos+1)]
        cons_vector = flatten(arrays)
        # Copy from -1 to boundary point at zero
        zero_eqn = self.create_boundary_equations(cons_vector, cons_vector, [(0, to_side)])
        kernel.add_equation(zero_eqn)
        output_eqns = self.create_boundary_equations(cons_vector, cons_vector, transfer_indices)
        kernel.add_equation(output_eqns)
        kernel.update_block_datasets(block)
        return kernel

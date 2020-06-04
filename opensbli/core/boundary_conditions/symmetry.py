from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from sympy import flatten, Matrix, sqrt


class SymmetryBC(BoundaryConditionBase):
    """ Applies a symmetry condition on the boundary, normal velocity components set/evaluate to zero.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'Symmetry'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        fd_metric = block.fd_metrics
        direction, side = self.direction, self.side
        direction_metric = Matrix(fd_metric[direction, :])
        normalisation = sqrt(sum([a**2 for a in direction_metric]))
        unit_normals = direction_metric/normalisation
        lhs_eqns = flatten(arrays)
        boundary_values = []
        rhs_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                contra_variant_vector = unit_normals.dot(Matrix(ar))
                transformed_vector = Matrix(ar).T - 2.*contra_variant_vector*unit_normals
                rhs_eqns += flatten(transformed_vector)
                # Later move this to an inviscid wall boundary condition
                transformed_vector = Matrix(ar).T - 1.*contra_variant_vector*unit_normals
                boundary_values += flatten(transformed_vector)
            else:
                rhs_eqns += [ar]
                boundary_values += [ar]
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        from_side_factor, to_side_factor = self.set_side_factor()

        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[direction][side]) + 1)]
        final_equations = self.create_boundary_equations(lhs_eqns, rhs_eqns, transfer_indices)
        transfer_indices = [tuple([0, to_side_factor])]
        final_equations += self.create_boundary_equations(lhs_eqns, boundary_values, transfer_indices)
        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)
        return kernel

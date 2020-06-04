from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from sympy import flatten
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.equation_types.opensbliequations import OpenSBLIEq


class ExtrapolationBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Extrapolation boundary condition. Pass order 0 for Zeroth order extrapolation
    and order=1 for linear extrapolation.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only.
    TODO Is it modify central"""

    def __init__(self, direction, side, order, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'Extrapolation'
        # Order of the extrapolation
        self.order = order
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        direction, side = self.direction, self.side
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        cons_vars = flatten(arrays)
        n_halos = abs(halos[self.direction][self.side])
        from_side_factor, to_side_factor = self.set_side_factor()

        if self.order == 0:  # Zeroth order extrapolation
            halo_points = [from_side_factor*i for i in range(0, n_halos+1)]
            for i in halo_points:
                equations = self.create_boundary_equations(cons_vars, cons_vars, [(i, to_side_factor)])
                kernel.add_equation(equations)
        elif self.order == 1:  # Linear extrapolation
            halo_points = [from_side_factor*i for i in range(1, n_halos+1)]
            equations = []
            if side == 0:
                for i in halo_points:
                    for term in cons_vars:
                        current, one, two = increment_dataset(term, direction, i), increment_dataset(term, direction, (i+1)*to_side_factor), increment_dataset(term, direction, (i+2)*to_side_factor)
                        equations += [OpenSBLIEq(current, 2.0*one - two)]
            elif side == 1:
                for i in halo_points:
                    for term in cons_vars:
                        current, one, two = increment_dataset(term, direction, i), increment_dataset(term, direction, (i-1)*from_side_factor), increment_dataset(term, direction, (i-2)*from_side_factor)
                        equations += [OpenSBLIEq(current, 2.0*one - two)]
            kernel.add_equation(equations)
        else:
            raise ValueError("Only zeroth order and linear extrapolation currently implemented.")
        kernel.update_block_datasets(block)
        return kernel

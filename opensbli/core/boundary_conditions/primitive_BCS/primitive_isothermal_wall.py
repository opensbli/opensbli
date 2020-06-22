from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.opensbliobjects import ConstantObject, DataObject
from sympy import Rational, S, Float, pprint
from opensbli.core.grid import GridVariable


class PrimitiveIsothermalWallBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Navier-Stokes specific boundary condition. Applies a no-slip viscous wall condition,
    velocity components are zero on the wall. Temperature is fixed with a prescribed wall temperature.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg Eq equations: Equation for conservative variable rhoE to set on the wall with constant temperature.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, wall_temperature, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'IsothermalWall'
        self.wall_temperature = wall_temperature
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        n_halos = abs(halos[self.direction][self.side])
        from_side_factor, to_side_factor = self.set_side_factor()
        velocity_components = [block.location_dataset('u%d' % i) for i in range(block.ndim)]
        # Set wall conditions, velocity components are zeroed
        wall_eqns = [OpenSBLIEq(x, Float(S.Zero)) for x in velocity_components]
        # Set wall temperature, density is left to be evaluated
        wall_eqns += [OpenSBLIEq(block.location_dataset('T'), self.wall_temperature)]
        kernel.add_equation(wall_eqns)
        kernel.update_block_datasets(block)
        return kernel

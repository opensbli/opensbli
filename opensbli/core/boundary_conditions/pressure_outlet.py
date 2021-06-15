from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.utilities.helperfunctions import dot
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.grid import GridVariable
from sympy import flatten
from opensbli.core.opensbliobjects import ConstantObject
from opensbli.core.kernel import ConstantsToDeclare


class PressureOutletBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Specified outlet (back) pressure boundary condition. Pressure is specified, density and velocity components
    are extrapolated from 1 point inside the boundary to the boundary point and the halos on that side.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg float back_pressure: Numerical value of back pressure to set on the outlet.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""
    def __init__(self, direction, side, back_pressure, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'PressureOutlet'
        self.back_pressure = back_pressure
        if side != 1:
            raise ValueError("Only implemented for side=1 outlet.")
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        NS = NSphysics(block)
        cons_vars = flatten([NS.density(), NS.momentum(), NS.total_energy()])
        rhs = flatten([NS.density(), NS.momentum()])
        n_halos = abs(halos[self.direction][self.side])
        local_rhoE = OpenSBLIEq(GridVariable('local_rhoE'), NS.total_energy())
        kernel.add_equation(local_rhoE)
        # from_side_factor, to_side_factor = self.set_side_factor()
        halo_points = [0] + [i for i in range(1, n_halos+1)]
        # Set rhoE based on the specified back pressure and extrapolated density and velocities
        bp = ConstantObject('back_pressure')
        ConstantsToDeclare.add_constant(bp)
        bp._value = self.back_pressure
        rhs += [bp/(NS.specific_heat_ratio()-1.0) + 0.5*dot(NS.momentum(), NS.momentum())/NS.density()]
        for i in halo_points:  # Extrapolate rho, rhou, rhov, rhow from one point inside the domain
            equations = self.create_boundary_equations(cons_vars, rhs, [(i, -1)])
            kernel.add_equation(equations)
        kernel.update_block_datasets(block)
        return kernel

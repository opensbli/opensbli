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
    velocity and temperature fluctuations are zero on the wall. Density fluctuations are left
    to change based on the continuity equation.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg Eq equations: Equation for conservative variable rhoE to set on the wall with constant temperature.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'IsothermalWall'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        n_halos = abs(halos[self.direction][self.side])
        from_side_factor, to_side_factor = self.set_side_factor()
        direction = self.direction
        velocity_components = [block.location_dataset('u%d' % i) for i in range(block.ndim)]
        # Set wall conditions, velocity components are zeroed
        wall_eqns = [OpenSBLIEq(x, Float(S.Zero)) for x in velocity_components]
        # Set wall temperature, density is left to be evaluated
        wall_eqns += [OpenSBLIEq(block.location_dataset('T'), 0.0)]
        kernel.add_equation(wall_eqns)
        # Set the halos
        halo_eqns = []
        # Reverse the velocity components over the wall
        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, n_halos + 1)]
        print(transfer_indices)
        for index, component in enumerate(velocity_components):
            for pair_index, pair in enumerate(transfer_indices):
                left = increment_dataset(component, direction, pair[0])
                right = increment_dataset(component, direction, pair[1])
                halo_eqns += [OpenSBLIEq(left, -right)]
        # Set the temperature in the halo points T[-i] = i*T[0]
        T = block.location_dataset('T')
        T1 = increment_dataset(T, direction, 1*to_side_factor)
        halo_eqns += [OpenSBLIEq(increment_dataset(T, direction, -i*to_side_factor), -i*T1) for i in range(1, n_halos+1)]
        # Wall pressure
        gama, Minf = ConstantObject('gama'), ConstantObject('Minf')
        rho, T = block.location_dataset('rho'), block.location_dataset('T')
        pw = OpenSBLIEq(GridVariable('pw'), rho*T/(gama*Minf**2))
        # Base flow wall pressure
        rhob, Tb = block.location_dataset('rhob'), block.location_dataset('Tb')
        pwb = OpenSBLIEq(GridVariable('pwb'), rhob*Tb/(gama*Minf**2))
        halo_eqns += [pw, pwb]
        # Set density in the halos
        for i in range(1, n_halos+1):
            T_halo, Tb_halo = increment_dataset(T, direction, -i*to_side_factor), increment_dataset(Tb, direction, i*to_side_factor)
            rhs = Minf**2 * gama * pwb.lhs/Tb_halo * (pw.lhs/pwb.lhs - T_halo/Tb_halo)
            rho_lhs = increment_dataset(rho, direction, -i*to_side_factor)
            halo_eqns += [OpenSBLIEq(rho_lhs, rhs)]
        kernel.add_equation(halo_eqns)
        print("Printing the wall boundary equations:\n")
        for eqn in kernel.equations:
            pprint(eqn)
        kernel.update_block_datasets(block)
        return kernel

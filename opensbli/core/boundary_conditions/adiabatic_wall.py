from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.opensbliobjects import ConstantObject
from sympy import Rational


class AdiabaticWallBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Adiabatic wall condition, zero gradient dT/dn = 0 over the boundary.
    # modified by Hiten Mulchandani (November 2019)

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'AdiabaticWall'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        NS = NSphysics(block)  # using Navier Stokes physics object, create conservative variables
        gama = NS.specific_heat_ratio()
        Minf = ConstantObject('Minf')
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        direction, side = self.direction, self.side
        n_halos = abs(halos[direction][side])
        from_side_factor, to_side_factor = self.set_side_factor()
        rho_wall = NS.density()
        halo_densities, momentum_wall, halo_momentum, halo_energy = [], [], [], []
        for index, momentum_comp in enumerate(NS.momentum()):
            velocity_wall = 0  # get u0, u1 components of velocity
            momentum_wall += [rho_wall*velocity_wall]  # get momentum at the wall
            halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, direction, 0), momentum_wall[index])]
        T_above = increment_dataset(NS.temperature(relation=True, conservative=True), direction, to_side_factor)
        T_above2 = increment_dataset(NS.temperature(relation=True, conservative=True), direction, 2*to_side_factor)
        T_above3 = increment_dataset(NS.temperature(relation=True, conservative=True), direction, 3*to_side_factor)
        # first-order approximation of dT/dy = 0
        # T_wall = T_above # first-order (to begin with)
        # second-order approximation of dT/dy = 0
        # T_wall = (4*T_above - T_above2)/3
        # Carpenter's one-sided fourth-order approach for dT/dy = 0
        T_wall = Rational(6, 11)*(3*T_above + Rational(1, 3)*T_above3 - Rational(3, 2)*T_above2)
        pressure_wall = rho_wall*T_wall/(gama*Minf*Minf)
        energy_wall = pressure_wall/(gama-1)
        halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), direction, 0), energy_wall)]
        for i in range(1, n_halos+1):
            extrapolated_density = increment_dataset(NS.density(), direction, to_side_factor*i)
            halo_densities += [OpenSBLIEq(increment_dataset(NS.density(), direction, from_side_factor*i), extrapolated_density)]
            momentum_above = []
            for index, momentum_comp in enumerate(NS.momentum()):
                momentum_above += [increment_dataset(momentum_comp, direction, to_side_factor*i)]
            for index, momentum_comp in enumerate(NS.momentum()):
                extrapolated_momentum = momentum_above[index]
                halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, direction, from_side_factor*i), -extrapolated_momentum)]
            extrapolated_energy = increment_dataset(NS.total_energy(), direction, to_side_factor*i)
            halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), direction, from_side_factor*i), extrapolated_energy)]
        kernel.add_equation(halo_densities)
        kernel.add_equation(halo_momentum)
        kernel.add_equation(halo_energy)
        kernel.update_block_datasets(block)
        return kernel

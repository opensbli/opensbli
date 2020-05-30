from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.grid import GridVariable
from sympy import Rational

class InletLawalBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Velocity is extrapolated from 1 point inside the boundary to the boundary point and the halos on that side. 
    All primitive and conservative variables at the inlet and in the halos are then constructed from the extrapolated
    velocity. At present, only zeroth order extrapolation has been implemented.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""
    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'InletLawal'
        if side != 0:
            raise ValueError("Only implemented for side=0 inlet.")
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        NS = NSphysics(block) # using Navier Stokes physics object, create conservative variables
        gama = NS.specific_heat_ratio()
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        direction, side = self.direction, self.side
        n_halos = abs(halos[direction][side])
        from_side_factor, to_side_factor = self.set_side_factor()
        local_rhoE = OpenSBLIEq(GridVariable('local_rhoE'), NS.total_energy())
        kernel.add_equation(local_rhoE)
        halo_velocities = []
        velocities = []
        # interpolate velocities from one point inside the domain
        for index, momentum_comp in enumerate(NS.momentum()):
            momentum, density = increment_dataset(momentum_comp, direction, 1), increment_dataset(NS.density(), direction, 1)
            halo_velocities += [OpenSBLIEq(GridVariable('u%d' % (index)), momentum/density)]
            velocity = GridVariable('u%d' % (index))
            velocities += [velocity**2]
        kernel.add_equation(halo_velocities)
        halo_densities = []
        halo_momentum = []
        halo_energy = []
        velocities_squared = sum(velocities)
        # set temperature in the halos based on interpolated velocities
        T_halo = 1.0-(Rational(1, 2)*(gama-1.0)*velocities_squared)
        # set density in the halos based on interpolated velocities and calculated temperature
        rho_halo = (1.0+Rational(1, 2)*(gama-1.0)*(velocities_squared/T_halo))**(-1.0/(gama-1.0))
        # set energy in the halos based on the calculated density, temperature and interpolated velocities
        energy_rhs = rho_halo*T_halo/(gama*(gama-1.0)) + Rational(1,2)*rho_halo*velocities_squared
        for i in range(0, n_halos+1):
            # set momentum in the halos based on the interpolated density and velocity
            for index, momentum_comp in enumerate(NS.momentum()):
                velocity = GridVariable('u%d' % (index))
                halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, direction, from_side_factor*i), rho_halo*velocity)]
            halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), direction, from_side_factor*i), energy_rhs)]
            halo_densities += [OpenSBLIEq(increment_dataset(NS.density(), direction, from_side_factor*i), rho_halo)]
        kernel.add_equation(halo_densities)
        kernel.add_equation(halo_momentum)
        kernel.add_equation(halo_energy)
        kernel.update_block_datasets(block)
        return kernel
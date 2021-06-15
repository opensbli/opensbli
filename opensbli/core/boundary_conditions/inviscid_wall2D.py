from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.utilities.helperfunctions import increment_dataset
from sympy.functions.elementary.piecewise import Piecewise
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.grid import GridVariable
from sympy import sqrt, atan2, cos, sin, Matrix, Rational


class InviscidWall2DBC(BoundaryConditionBase):  # A.A. Lawal's method, see pg. 88 of his PhD thesis.
    """ Applies a slip condition on the boundary; normal velocity components evaluate to zero.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.bc_name = 'InviscidWall2D'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        NS = NSphysics(block)
        gama = NS.specific_heat_ratio()
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        direction, side = self.direction, self.side
        tangent_direction = 1-direction
        n_halos = abs(halos[direction][side])
        from_side_factor, to_side_factor = self.set_side_factor()
        dimension = block.ndim
        # get x and y at the wall for the tangent angle calculation
        x0_wall = block.location_dataset('x0')
        x1_wall = block.location_dataset('x1')
        # get x and y at one grid point away from the wall for the grid size calculation
        x0 = increment_dataset(block.location_dataset('x0'), direction, to_side_factor)
        x1 = increment_dataset(block.location_dataset('x1'), direction, to_side_factor)
        # get density, momentum, and pressure at one grid point from the wall
        density_above = increment_dataset(NS.density(relation=True), direction, to_side_factor)
        momentum_above = []
        for index, momentum_comp in enumerate(NS.momentum()):
            momentum_above += [increment_dataset(momentum_comp, direction, to_side_factor)]
        pressure_above = increment_dataset(NS.pressure(relation=True, conservative=True), direction, to_side_factor)
        energy_above = increment_dataset(NS.total_energy(), direction, to_side_factor)
        # to obtain flow-field properties in the tangential direction
        tangent_indices = [-1, 0, 1]
        # to calculate the gradient of the geometry using the central difference (CD) scheme
        gradient_indices = [-2] + tangent_indices + [2]
        # get x and y at the tangent indices
        x, y, x_wall, y_wall = [], [], [], []
        for index, gradient_index in enumerate(gradient_indices):
            x += [increment_dataset(x0, tangent_direction, gradient_index)]
            y += [increment_dataset(x1, tangent_direction, gradient_index)]
            x_wall += [increment_dataset(x0_wall, tangent_direction, gradient_index)]
            y_wall += [increment_dataset(x1_wall, tangent_direction, gradient_index)]
        # to calculate the gradient of the wall using the second order CD scheme
        delta_y, delta_x = [], []
        for index, tangent_index in enumerate(tangent_indices):
            index_1 = tangent_index+1
            index_2 = index_1+2
            delta_y += [y_wall[index_2]-y_wall[index_1]]
            delta_x += [x_wall[index_2]-x_wall[index_1]]
        # calculate density, velocity, and pressure components at the (i-1), i, and (i+1) grid points
        theta_calc, density_calc, rho, velocity_calc, pressure_calc, p = [], [], [], [], [], []
        for index_1, tangent_index in enumerate(tangent_indices):
            y_arg = delta_y[index_1]
            x_arg = delta_x[index_1]
            theta_calc += [OpenSBLIEq(GridVariable('theta_%d' % (1+tangent_index)), atan2(y_arg, x_arg))]  # tangent angle at the wall
            density_calc += [OpenSBLIEq(GridVariable('rho_%d' % (1+tangent_index)), increment_dataset(density_above, tangent_direction, tangent_index))]
            rho += [GridVariable('rho_%d' % (1+tangent_index))]
            for index_2, momentum_comp in enumerate(momentum_above):
                momentum_calc = increment_dataset(momentum_comp, tangent_direction, tangent_index)
                velocity_calc += [OpenSBLIEq(GridVariable('u_%d%d' % (index_2, 1+tangent_index)), momentum_calc/rho[index_1])]
            pressure_calc += [OpenSBLIEq(GridVariable('p_%d' % (1+tangent_index)), increment_dataset(pressure_above, tangent_direction, tangent_index))]
            p += [GridVariable('p_%d' % (1+tangent_index))]
        kernel.add_equation(theta_calc)
        kernel.add_equation(density_calc)
        kernel.add_equation(velocity_calc)
        kernel.add_equation(pressure_calc)
        # calculate tangential velocity at the (i-1), i, and (i+1) grid points
        theta, tangent_velocity, ut = [], [], []
        unit_normals = [[0]*dimension for i in range(len(tangent_indices))]
        for index_1, tangent_index in enumerate(tangent_indices):
            theta += [GridVariable('theta_%d' % (1+tangent_index))]
            unit_normals[index_1] = [cos(theta[index_1]), sin(theta[index_1])]
            velocity = []
            for index_2, momentum_comp in enumerate(momentum_above):
                velocity += [GridVariable('u_%d%d' % (index_2, 1+tangent_index))]
            unit_normals_array = Matrix(unit_normals[index_1]).T
            tangent_velocity += [OpenSBLIEq(GridVariable('ut_%d' % (1+tangent_index)), unit_normals_array.dot(Matrix(velocity)))]
            ut += [GridVariable('ut_%d' % (1+tangent_index))]
        kernel.add_equation(tangent_velocity)
        # calculate the size of the wall grid point
        grid_size = Piecewise((sqrt((y[2]-y[1])**2+(x[2]-x[1])**2), theta[1] >= 0), (sqrt((y[3]-y[2])**2+(x[3]-x[2])**2), True))
        # calculate the difference in density, tangential velocity, and pressure at the ith grid point
        density_delta = Piecewise(((rho[1]-rho[0]), theta[1] >= 0), ((rho[1]-rho[2]), True))
        tangent_velocity_delta = Piecewise(((ut[1]-ut[0]), theta[1] >= 0), ((ut[1]-ut[2]), True))
        pressure_delta = Piecewise(((p[1]-p[0]), theta[1] >= 0), ((p[1]-p[2]), True))
        # calculate the location of the grid point normal to the wall
        normal_grid_point = (x1-x1_wall)*sin(theta[1])/grid_size
        # calculate the density, tangent velocity, and pressure normal to ith grid point
        density_wall = rho[1]-(normal_grid_point*density_delta)
        tangent_velocity_prime = ut[1]-(normal_grid_point*tangent_velocity_delta)
        pressure_wall = p[1]-(normal_grid_point*pressure_delta)
        halo_densities, momentum_wall, halo_momentum, halo_energy = [], [], [], []
        # assign density, momentum, and energy at the wall
        halo_densities += [OpenSBLIEq(increment_dataset(NS.density(), direction, 0), density_wall)]
        for index, momentum_comp in enumerate(NS.momentum()):
            projection = unit_normals[1][index]  # get tangent angle at the wall
            velocity_wall = tangent_velocity_prime*projection  # get u0, u1 components of velocity
            momentum_wall += [density_wall*velocity_wall]  # get momentum at the wall
            halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, direction, 0), momentum_wall[index])]
        energy_wall = pressure_wall/(gama-1) + Rational(1, 2)*density_wall*(tangent_velocity_prime**2)
        halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), direction, 0), energy_wall)]
        # assign density, momentum, and energy in the halos (see section 8.2.1 in CFD textbook by J. Blazek for details)
        for i in range(1, n_halos+1):
            extrapolated_density = (i+1)*density_wall - i*density_above
            halo_densities += [OpenSBLIEq(increment_dataset(NS.density(), direction, from_side_factor*i), extrapolated_density)]
            for index, momentum_comp in enumerate(NS.momentum()):
                extrapolated_momentum = (i+1)*momentum_wall[index] - i*momentum_above[index]
                halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, direction, from_side_factor*i), extrapolated_momentum)]
            extrapolated_energy = (i+1)*energy_wall - i*energy_above
            halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), direction, from_side_factor*i), extrapolated_energy)]
        kernel.add_equation(halo_densities)
        kernel.add_equation(halo_momentum)
        kernel.add_equation(halo_energy)
        kernel.update_block_datasets(block)
        return kernel

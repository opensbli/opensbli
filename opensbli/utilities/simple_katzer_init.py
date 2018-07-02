from sympy import Piecewise
import numpy as np
from opensbli.initialisation import GridBasedInitialisation
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.core.grid import GridVariable
from opensbli.core.kernel import Kernel
from opensbli.equation_types.opensbliequations import OpenSBLIEq


class Initialise_Katzer(GridBasedInitialisation):
    """ Generates the initialiastion equations for the boundary-layer profile.
    :arg list bl_directions: Integer values of the problem directions.
    arg list coordinate_evaluations: User specified formulae to generate a grid."""
    def __new__(cls, bl_directions, coordinate_evaluations=None):
        ret = super(Initialise_Katzer, cls).__new__(cls)
        ret.coordinates = [x[1] for x in bl_directions]
        ret.bl_directions = bl_directions
        ret.coordinate_evaluations = coordinate_evaluations
        ret.equations = []
        return ret

    def check_inputs(self, block):
        bl_directions = [x[0] for x in self.bl_directions]
        if sum(bl_directions) < 1:
            raise ValueError("Provide the directions to apply a boundary layer profile in.")
        if len(bl_directions) != block.ndim:
            raise ValueError("The list of polynomial directions must match the dimensions of the problem.")
        return

    def spatial_discretisation(self, block):
        self.equations = []
        self.block = block
        self.idxs = block.grid_indexes
        self.check_inputs(block)
        # Check if user has passed equations to evaluate coordinates, and add them to the kernel
        if self.coordinate_evaluations:
            self.equations += self.coordinate_evaluations
        self.initial = self.generate_initial_condition()
        # Add polynomial equations to initialise the solution
        self.equations += self.eqns
        self.equations = block.dataobjects_to_datasets_on_block(self.equations)
        kernel1 = Kernel(block, computation_name="Grid_based_initialisation%d" % self.order)
        kernel1.set_grid_range(block)
        schemes = block.discretisation_schemes
        for d in range(block.ndim):
            for sc in schemes:
                if schemes[sc].schemetype == "Spatial":
                    kernel1.set_halo_range(d, 0, schemes[sc].halotype)
                    kernel1.set_halo_range(d, 1, schemes[sc].halotype)
        kernel1.add_equation(self.equations)
        kernel1.update_block_datasets(block)
        self.Kernels = [kernel1]
        return

    def generate_initial_condition(self):
        bl_directions = [x[0] for x in self.bl_directions]
        # Load pre-computed profiles for rhou, rhov, and T
        rhou_coeffs, rhov_coeffs, T_coeffs = self.load_similarity()
        # x1 coordinate at the edge of the boundary-layer near the inlet
        bl_edge = 2.316

        if sum(bl_directions) == 1:  # 2D Katzer and 3D spanwise periodic Katzer, boundary layer in one direction
            dire = [i for i, x in enumerate(self.bl_directions) if x[0]][0]
            self.generate_one_wall_equations([rhou_coeffs, rhov_coeffs, T_coeffs], dire, bl_edge)
        else:
            raise NotImplementedError("Only one direction can have a boundary-layer in this initialistion.")
        return

    def form_equation(self, name, coefficients, direction, edge, free_stream):
        """ Creates the piecewise equations for the cases of 2D and 3D span-periodic boundary-layer profiles.

        :arg string name: Name of the variable.
        :arg ndarray coefficients: Coefficients for the polynomial fit.
        :arg int direction: Spatial direction to apply the equation to.
        :arg int edge: Grid index for the edge of the boundary-layer.
        :arg float free_stream: Free stream value of this variable.
        returns: Eq: eqn: OpenSBLI equation to add to the initialisation kernel."""
        powers = [i for i in range(np.size(coefficients))][::-1]
        eqn = sum([coeff*self.coordinates[direction]**power for (coeff, power) in zip(coefficients, powers)])
        eqn = OpenSBLIEq(GridVariable('%s' % name), Piecewise((eqn, self.coordinates[direction] < edge), (free_stream, True)))
        return eqn

    def generate_one_wall_equations(self, coeffs, direction, edge):
        """ Generates the equations for 2D SBLI and 3D span-periodic cases.

        :arg list data: Profile arrays for rhou0, rhou1 and temperature.
        :arg list coeffs: Coefficients for the polynomial fits.
        :arg list direction: Direction normal to the wall.
        :arg int edge: Coordinate for the edge of the boundary-layer."""
        self.eqns = []
        rhou0_eqn = self.form_equation('rhou0', coeffs[0], direction, edge, 1.0)
        rhou1_eqn = self.form_equation('rhou1', coeffs[1], direction, edge, 0.005650016)
        T_eqn = self.form_equation('T', coeffs[2], direction, edge, 1.0)
        # Set conservative values
        rho, rhou0, rhou1, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('T')
        rho_eqn = OpenSBLIEq(rho, 1.0/T)
        rho_store = OpenSBLIEq(DataObject('rho'), rho)
        rhou0_store = OpenSBLIEq(DataObject('rhou0'), rhou0)
        rhou1_store = OpenSBLIEq(DataObject('rhou1'), rhou1)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = OpenSBLIEq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2)/rho)
        self.eqns += [rhou0_eqn, rhou1_eqn, T_eqn, rho_eqn, rho_store, rhou0_store, rhou1_store, rhoE_store]
        if self.block.ndim == 3:  # Periodic case, rhow = 0
            self.eqns += [OpenSBLIEq(DataObject('rhou2'), 0.0)]
        return

    def load_similarity(self):
        """ Pre-computed coefficients for a similarity solution to the compressible
        boundary-layer eqautions."""
        rhou = np.array([-0.000661, 0.456679, -0.378541, 1.829143, -4.060705, 5.261448, -4.002925, 1.796903, -0.469165, 0.066022, -0.003878])[::-1]
        rhov = np.array([-0.000019, 0.000759, -0.006145, 0.029622, -0.061263, 0.072858, -0.051201, 0.021523, -0.005326, 0.000717, -0.000041])[::-1]
        T = np.array([1.675785, 0.011027, -0.346234, 0.031632, 0.310312, -0.838071, 0.945210, -0.530435, 0.158965, -0.024521, 0.001537])[::-1]
        return rhou, rhov, T

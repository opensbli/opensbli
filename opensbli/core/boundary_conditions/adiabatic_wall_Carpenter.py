from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from sympy import sqrt, Matrix, pprint
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.opensbliobjects import ConstantObject
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.core.grid import GridVariable


class AdiabaticWall_CarpenterBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Adiabatic wall condition, zero gradient dT/dn = 0 over the boundary. (G.N. Coleman)
    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'AdiabaticWall_Carpenter'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme

        self.deriv_indices = [[i for i in range(6)], [-i for i in range(6)]]
        self.deriv_coeffs = [[self.carp_coefficients()[0, i] for i in range(6)], [-self.carp_coefficients()[0, i] for i in range(6)]]
        return

    def carp_coefficients(self):
        """ Computes the finite-difference coefficients for the 1st order one-sided Carpenter wall boundary derivative.
        :returns: Matrix: bc4: Matrix of stencil coefficients."""
        R1 = -(2177.0*sqrt(295369.0)-1166427.0)/25488.0
        R2 = (66195.0*sqrt(53.0)*sqrt(5573.0)-35909375.0)/101952.0

        al4_0 = [-(216.0*R2+2160.0*R1-2125.0)/12960.0, (81.0*R2+675.0*R1+415.0)/540.0, -(72.0*R2+720.0*R1+445.0)/1440.0, -(108.0*R2+756.0*R1+421.0)/1296.0]
        al4_1 = [(81.0*R2+675.0*R1+415.0)/540.0, -(4104.0*R2+32400.0*R1+11225.0)/4320.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(216.0*R2+2160.0*R1+655.0)/4320.0]
        al4_2 = [-(72.0*R2+720.0*R1+445.0)/1440.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(4104.0*R2+32400.0*R1+12785.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0]
        al4_3 = [-(108.0*R2+756.0*R1+421.0)/1296.0, -(216.0*R2+2160.0*R1+655.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0, -(216.0*R2+2160.0*R1-12085.0)/12960.0]
        al4 = Matrix([al4_0, al4_1, al4_2, al4_3])

        ar4_0 = [(-1.0)/2.0, -(864.0*R2+6480.0*R1+305.0)/4320.0, (216.0*R2+1620.0*R1+725.0)/540.0, -(864.0*R2+6480.0*R1+3335.0)/4320.0, 0.0, 0.0]
        ar4_1 = [(864.0*R2+6480.0*R1+305.0)/4320.0, 0.0, -(864.0*R2+6480.0*R1+2315.0)/1440.0, (108.0*R2+810.0*R1+415.0)/270.0, 0.0, 0.0]
        ar4_2 = [-(216.0*R2+1620.0*R1+725.0)/540.0, (864.0*R2+6480.0*R1+2315.0)/1440.0, 0.0, -(864.0*R2+6480.0*R1+785.0)/4320.0, -1.0/12.0, 0.0]
        ar4_3 = [(864.0*R2+6480.0*R1+3335.0)/4320.0, -(108.0*R2+810.0*R1+415.0)/270.0, (864.0*R2+6480.0*R1+785.0)/4320.0, 0.0, 8.0/12.0, -1.0/12.0]
        ar4 = Matrix([ar4_0, ar4_1, ar4_2, ar4_3])
        # Form inverse and convert to rational
        al4_inv = al4.inv()
        bc4 = al4_inv*ar4
        return bc4

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        wall_eqns = []
        direction, side = self.direction, self.side
        from_side_factor, to_side_factor = self.set_side_factor()
        NS = NSphysics(block)
        for ar in arrays:
            if isinstance(ar, list):  # Set velocity components to zero on the wall
                rhs = [0 for i in range(len(ar))]
                wall_eqns += [OpenSBLIEq(x, y) for (x, y) in zip(ar, rhs)]
        # Calculate Temperature at the wall associated with  dT/dy = 0
        coeffs, grid_indices = self.deriv_coeffs[side], self.deriv_indices[side]
        # Evaluate temperature values at the wall and 5 points above/below it
        temperature_values = [GridVariable('T%d' % i) for i in range(len(grid_indices))]
        temperature_evaluations = [OpenSBLIEq(temperature_values[i], increment_dataset(NS.temperature(relation=True, conservative=True), self.direction, i*to_side_factor))
                                   for i in range(1, len(temperature_values))]

        wall_eqns += temperature_evaluations
        deriv_equation = sum([c*temperature_values[i+1] for i, c in enumerate(coeffs[1:])])
        # Rearrange the equation to solve for the wall temperature
        deriv_equation = OpenSBLIEq(temperature_values[0], -1.0*deriv_equation/(coeffs[0]))  # set wall value of dT/dy to 0...
        wall_eqns += [deriv_equation]

        # Set energy on the wall based on this calculated wall temperature and the density that comes from the continuity equation
        Minf, gama = ConstantObject('Minf'), ConstantObject('gama')
        wall_eqns += [OpenSBLIEq(NS.total_energy(), NS.density()*temperature_values[0] / (Minf*Minf*gama*(gama-1.0)))]
        kernel.add_equation(wall_eqns)
        # Print out the current equations for the adiabatic wall conditio
        print("Printing equations for Adiabatic wall boundary condition (Carp-4):")  # <<<
        for eqn in kernel.equations:
            pprint(eqn)
        # exit()
        kernel.update_block_datasets(block)
        return kernel

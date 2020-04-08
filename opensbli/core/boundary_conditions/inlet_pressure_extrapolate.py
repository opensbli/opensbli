from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase, ModifyCentralDerivative
from opensbli.core.boundary_conditions.Carpenter_scheme import Carpenter
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.core.grid import GridVariable
from sympy import Abs, flatten, GreaterThan, Piecewise
from sympy.functions.elementary.piecewise import ExprCondPair


class InletPressureExtrapolateBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Navier-Stoke boundary. Inlet boundary condition, local velocity normal to the boundary
    (for cartesian grid only) is compared to the local speed of sound,
    pressure is extrapolated out into the halos if subsonic.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'InletPressureExtrapolate'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        if side != 0:
            raise ValueError("Only implemented this BC for inlet side 0.")
        return

    def apply(self, arrays, block):
        direction = self.direction
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        # Using Navier Stokes physics object, create conservative variables
        NS = NSphysics(block)
        cons_vars = [NS.density(), NS.momentum(), NS.total_energy()]
        # Evaluation of pressure, density, speed of sound on the boundary
        pb, rhob, ab = GridVariable('pb'), GridVariable('rhob'), GridVariable('ab')
        gama = NS.specific_heat_ratio()
        grid_vels = [GridVariable('ub%d' % i) for i in range(block.ndim)]
        grid_vels_sq = [i**2 for i in grid_vels]
        eqns = [OpenSBLIEq(rhob, NS.density())]
        eqns += [OpenSBLIEq(grid_vels[i], Abs(u/NS.density())) for i, u in enumerate(NS.momentum())]
        eqns += [OpenSBLIEq(pb, (gama-1)*(flatten(arrays)[-1] - 0.5*rhob*sum(flatten(grid_vels_sq))))]
        eqns += [OpenSBLIEq(ab, (gama*pb/rhob)**0.5)]
        kernel.add_equation(eqns)
        locations = [-1, 0]
        inlet_vel = grid_vels[direction]
        # Conditions to be set at the boundary
        for lhs in flatten(cons_vars):
            ecs = []
            rhs_values = [increment_dataset(lhs, direction, value) for value in locations]
            ecs += [ExprCondPair(rhs_values[0], GreaterThan(inlet_vel, ab))]
            ecs += [ExprCondPair(rhs_values[1], True)]
            kernel.add_equation(OpenSBLIEq(lhs, Piecewise(*ecs, **{'evaluate': False})))
        # Conditions set in the halos in rhoE
        locations = [-i-1 for i in range(abs(halos[0][0]))]
        lhs_rhoE = [increment_dataset(NS.total_energy(), direction, value) for value in locations]
        for i, lhs in enumerate(lhs_rhoE):
            ecs = []
            ecs += [ExprCondPair(lhs, GreaterThan(inlet_vel, ab))]  # lhs == rhs
            ecs += [ExprCondPair(NS.total_energy(), True)]
            kernel.add_equation(OpenSBLIEq(lhs, Piecewise(*ecs, **{'evaluate': False})))
        kernel.update_block_datasets(block)
        return kernel

#!/usr/bin/env python
from opensbli import *
import copy
from opensbli.utilities.katzer_init import Initialise_Katzer
from opensbli.utilities.helperfunctions import substitute_simulation_parameters, output_hdf5, print_iteration_ops
from opensbli.utilities.oblique_shock import ShockConditions
from sympy import tan, pi, tanh, sinh, cosh, exp, cos, sin
import time_averaging
# settings to turn on statistics gathering & read from restart file
stats = True
restart = True
stats_class = []
if stats:
    stats_class = time_averaging.get_stats_classes()

# Declare constant values
restart_iteration_no = symbols("restart_iteration_no", **{'cls': ConstantObject})
restart_iteration_no.datatype = Int()
CTD.add_constant([restart_iteration_no])

ndim = 3
sc1 = "**{\'scheme\':\'Teno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
# Substitutions
substitutions = [stress_tensor, heat_flux]
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Instatiate equation classes
eq = EinsteinEquation()
base_eqns = [mass, momentum, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity]
# Expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
# Expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)
# Local dictionary for the block
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Metric transformation
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 2)

# Create the Ducros equations for the shock sensor
SS = ShockSensor()
shock_sensor, sensor_array = SS.ducros_equations(block, coordinate_symbol, metriceq)
store_sensor = True
teno_order = 6
Avg = RoeAverage([0, 1])
LLF = LLFTeno(teno_order, formulation='adaptive', averaging=Avg, sensor=sensor_array, store_sensor=store_sensor)
schemes = {}
schemes[LLF.name] = LLF
fns = 'u0 u1 u2'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
block.set_discretisation_schemes(schemes)

# Create the body force after specifying a time stepping scheme
A, B, xF, yF, bta, omega, dt = symbols('A B xF yF bta omega dt', **{'cls': ConstantObject})
forcing_const = ["A", "B", "xF", "yF", "bta", "omega", "dt"]

current_iter = block.get_temporal_schemes[0].iteration_number
x0, x1, x2 = symbols('x0 x1 x2', **{'cls': DataObject})
# damping = (1 - exp(-x1 / 0.5))
body_force = Eq(DataObject('BF'), A*exp(-((x0-xF)**2+(x1-yF)**2))*cos(bta*x2)*sin(omega*dt*current_iter))

# Add forcing as an acoustic term to the continuity equation
for i, eq in enumerate(base_eqns):
    if i == 0:
        base_eqns[i] = Eq(eq.lhs, eq.rhs + body_force.lhs)

# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()
for eqn in base_eqns:
    simulation_eq.add_equations(eqn)
for eqn in constituent_eqns:
    constituent.add_equations(eqn)

# Apply metric transformation to the simulation equations
simulation_eq.apply_metrics(metriceq)


# Define a shock sensor for the TENO schemee
# Add the shock sensor to constituent relations
constituent.add_equations(shock_sensor)
# Add the body forcing term to the constituent relations
constituent.add_equations(body_force)

# Boundary conditions onwards
boundaries = [[0, 0] for t in range(ndim)]
# Left pressure extrapolation at x= 0, inlet conditions
direction, side = 0, 0
boundaries[direction][side] = InletPressureExtrapolateBC(direction, side, scheme=ReducedAccess())
# Right extrapolation at outlet
direction, side = 0, 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0, scheme=ReducedAccess())
# Bottom no-slip isothermal wall
local_dict['Lx1'] = ConstantObject('Lx1')
local_dict['by'] = ConstantObject('by')
direction, side = 1, 0
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns, scheme=ReducedAccess())
# Top dirichlet shock generator condition
direction, side = 1, 1
wave_angle = 44.6607551
xmach = 1.5
gamma = 1.4
pre_shock = (1.0, 1.0, 0.00466654053208844, (1.0/(gamma*xmach**2))/(gamma-1.0) + 0.5*(1.0*1.0**2 + 0.00466654053208844**2))
OS = ShockConditions(wave_angle, xmach, gamma)
post_shock = OS.conservative_post_shock_conditions(1.0)

rho = parse_expr("Eq(DataObject(rho), Piecewise((%.15f, DataObject(x0)>20.0), (%.15f, True)))" % (post_shock[0], pre_shock[0]), local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), Piecewise((%.15f, DataObject(x0)>20.0), (%.15f, True)))" % (post_shock[1], pre_shock[1]), local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), Piecewise((%.15f, DataObject(x0)>20.0), (%.15f, True)))" % (post_shock[2], pre_shock[2]), local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), Piecewise((%.15f + 0.5*DataObject(rhou2)**2 / DataObject(rho), DataObject(x0)>20.0), (%.15f + 0.5*DataObject(rhou2)**2 / DataObject(rho), True)))" % (post_shock[3], pre_shock[3]), local_dict=local_dict)
upper_eqns = [rho, rhou0, rhou1, rhoE]

boundaries[direction][side] = DirichletBC(direction, side, upper_eqns, scheme=ReducedAccess())
# Periodic direction 2
direction = 2
for side in [0,1]:
    boundaries[direction][side] = PeriodicBC(direction, side)
block.set_block_boundaries(boundaries)

# Perform initial condition
Re, xMach, Tinf = 750.0, 1.5, 202.17
# Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
grid_const = ["Lx1", "by"]
for con in grid_const:
    local_dict[con] = ConstantObject(con)
gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
gridx2 = parse_expr("Eq(DataObject(x2), block.deltas[2]*block.grid_indexes[2])", local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1, gridx2]
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1')), (False, DataObject('x2'))]
n_poly_coefficients = 50

initial = Initialise_Katzer(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_evaluations=coordinate_evaluation)
kwargs = {'iotype': "Write"}
h5 = iohdf5(arrays=simulation_eq.time_advance_arrays, save_every=10000, **kwargs)
h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('x2'), DataObject('D11')])
if store_sensor:
    h5.add_arrays([DataObject('TENO')])
block.setio(h5)

if stats:
    stats_arrays = time_averaging.get_arrays()
    kwargs = {'iotype': "Write", "name": "stats_output.h5"}
    h5_stats = iohdf5(**kwargs)
    h5_stats.add_arrays(stats_arrays)
    block.setio(h5_stats)

# Set equations on the block and discretise
block.set_equations([simulation_eq, constituent, initial, metriceq] + stats_class)
block.discretise()
# Create an algorithm and write the OPS C code
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
# Substitute the simulation parameters
constants = ['gama', 'Minf', 'Pr', 'Re', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2',
             'Delta0block0', 'Delta1block0', 'Delta2block0', 'SuthT', 'RefT', 'eps', 'Lx1', 'by', 'A', 'bta', 'omega', 'xF', 'yF', 'teno_a1', 'teno_a2', 'epsilon', 'restart_iteration_no']
values = ['1.4', '1.5', '0.72', '750.0', '1.3809973268575328', '0.025', '300000', '500', '200', '100',
          '375.0/(block0np0-1)', '140.0/(block0np1-1)', '27.32/(block0np2)', '110.4', '202.17', '1e-30', '140.0', '5.0', '2.5e-3', '0.23', '0.1011', '20.0', '4.0', '9.5', '3.5', '1.0e-16', '0']
substitute_simulation_parameters(constants, values)
print_iteration_ops()

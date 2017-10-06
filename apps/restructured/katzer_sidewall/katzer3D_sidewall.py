#!/usr/bin/env python
from opensbli.core import *
from opensbli.utilities.katzer_init import *
from opensbli.core.weno_opensbli import SimpleAverage
from opensbli.core.teno import LLFTeno
from opensbli.physical_models.euler_eigensystem import *
from opensbli.initialisation import GridBasedInitialisation
import copy

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

block = SimulationBlock(ndim, block_number=0)

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

local_dict['Lx1'], local_dict['Lx2'] = ConstantObject('Lx1'), ConstantObject('Lx2')
local_dict['by'], local_dict['bz'] = ConstantObject('by'), ConstantObject('bz')

gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
gridx2 = parse_expr("Eq(DataObject(x2), Lx2*sinh(bz*block.deltas[2]*block.grid_indexes[2]/Lx2)/sinh(bz))", local_dict=local_dict)

x_loc = parse_expr("Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)


rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5* d *(u0**2+u1**2))", local_dict=local_dict)

boundaries = [[0, 0] for t in range(ndim)]
# Left pressure extrapolation at x= 0, inlet conditions
direction = 0
side = 0
boundaries[direction][side] = InletTransferBoundaryConditionBlock(direction, side)
# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = OutletTransferBoundaryConditionBlock(direction, side)
# Bottom no-slip isothermal wall
direction = 1
side = 0
# boundaries[direction][side] = SymmetryBoundaryConditionBlock(direction, side)
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)
# Isothermal wall condition
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBoundaryConditionBlock(direction, side, wall_eqns, local_dict)
# Top dirichlet shock generator condition
direction = 1
side = 1
rho = parse_expr("Eq(DataObject(rho), Piecewise((1.129734572, (x0)>40.0), (1.00000596004, True)))", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), Piecewise((1.0921171, (x0)>40.0), (1.00000268202, True)))", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), Piecewise((-0.058866065, (x0)>40.0), (0.00565001630205, True)))", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), Piecewise((1.0590824, (x0)>40.0), (0.94644428042, True)))", local_dict=local_dict)

upper_eqns = [x_loc, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBoundaryConditionBlock(direction, side, upper_eqns)

direction = 2
side = 0
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBoundaryConditionBlock(direction, side, wall_eqns, local_dict)

direction = 2
side = 1
boundaries[direction][side] = SymmetryBoundaryConditionBlock(direction, side)

# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

for eqn in base_eqns:
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)
teno_order = 5
Euler_eq = EulerEquations(ndim)
ev_dict, LEV_dict, REV_dict = Euler_eq.generate_eig_system()
Avg = SimpleAverage([0, 1])
LLF = LLFTeno(ev_dict, LEV_dict, REV_dict, teno_order, ndim, Avg)
schemes = {}
schemes[LLF.name] = LLF
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

metriceq = MetricsEquation()
# Stretched in wall normal and spanwise directions
metriceq.genreate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (True, False)], 2)
simulation_eq.apply_metrics(metriceq)

block.set_block_boundaries(boundaries)

sim_eq = copy.deepcopy(simulation_eq)
CR = copy.deepcopy(constituent)

# Perform initial condition
coordinate_evaluation = [gridx0, gridx1, gridx2]
# Call the new polynomial based katzer initialisation, stretch factor 3 with 17 coefficients for the polynomial
init_katzer = Initialise_Katzer([600, 500, 150], [400.0, 115.0, 50.0], [1, 2], [3.0, 3.0], 17, block, coordinate_evaluation)
initial = init_katzer.initial

# Set equations on the block and discretise
block.set_equations([sim_eq, CR, initial, metriceq])
block.set_discretisation_schemes(schemes)
block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

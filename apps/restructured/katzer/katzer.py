#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *
from sympy import *
from opensbli.initialisation import *
from opensbli.utilities.gridgen import *

ndim = 2
weno_order = '5Z'
# weno_order = 3
Euler_eq = EulerEquations(ndim)
ev_dict, LEV_dict, REV_dict = Euler_eq.generate_eig_system()
Avg = SimpleAverage([0, 1])
LLF = LLFCharacteristic(ev_dict, LEV_dict, REV_dict, weno_order, ndim, Avg)

sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))"% sc1
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1

stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))" 
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
# Substitutions
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr","gama", "Minf"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# metriceq = MetricsEquation()
# metriceq.genreate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)



# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
#mu_B0[OPS_ACC1(0,0)] = pow(T_B0[OPS_ACC0(0,0)], 1.5)*(1.0+SutherlandT/ReferenceT)/(T_B0[OPS_ACC0(0,0)]+SutherlandT/ReferenceT);

viscosity = "Eq(mu, (T**(1.5)*(1.0+110.4/288.0)/(T+110.4/288.0)))"


simulation_eq = SimulationEquations()
eq = EinsteinEquation()
eqns = eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# simulation_eq.apply_metrics(metriceq)

constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(speed_of_sound, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

constituent.add_equations(EinsteinEquation().expand(temperature, ndim, coordinate_symbol, substitutions, constants))

eqns = eq.expand(viscosity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

schemes = {}
cent = Central(4)
schemes[cent.name] = cent
schemes[LLF.name] = LLF
rk = RungeKutta(3)
schemes[rk.name] = rk

block = SimulationBlock(ndim, block_number = 0)

# block.sbli_rhs_discretisation = True

local_dict = {"block" : block, "GridVariable" : GridVariable, "DataObject" : DataObject}

constants += ['Lx1', 'stretchfactor']
for con in constants:
    local_dict[con] = ConstantObject(con)

# Initial conditions
# x0 = parse_expr("Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
# x1 = parse_expr("Eq(GridVariable(x1), block.deltas[1]*block.grid_indexes[1])", local_dict=local_dict)

gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), block.deltas[1]*block.grid_indexes[1])", local_dict=local_dict)

# gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(stretchfactor*block.grid_indexes[1]/block.ranges[1][1])/sinh(stretchfactor))", local_dict=local_dict)

const_dict = {ConstantObject('stretchfactor') : 1.0, ConstantObject('Lx1') : 10, block.grid_indexes[1] : 2, block.ranges[1][1] : 10}


rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)",local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5* d *(u0**2+u1**2))", local_dict=local_dict)



d_in = parse_expr("Eq(GridVariable(d), 1.0)", local_dict=local_dict)
u0_in = parse_expr("Eq(GridVariable(u0), 1.0)", local_dict=local_dict)
u1_in = parse_expr("Eq(GridVariable(u1), 0.0)", local_dict=local_dict)
p_in = parse_expr("Eq(GridVariable(p), 1.0/(gama*Minf**2.0))", local_dict=local_dict) # p = 1/(gama*Minf^2)

initial = GridBasedInitialisation()
initial.add_equations([ gridx0, gridx1])

boundaries = [[0,0] for t in range(ndim)]
# Left Dirchlet at x= 0, inlet conditions
### Inflow conditions
direction = 0
side = 0
# inlet_eq = [d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE]
# boundaries[direction][side] = DirichletBoundaryConditionBlock(direction, side, inlet_eq)
boundaries[direction][side] = InletPressureExtrapolateBoundaryConditionBlock(direction, side)

# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = LinearExtrapolateBoundaryConditionBlock(direction, side)

# Bottom inviscid wall
direction = 1; side = 0 
# boundaries[direction][side] = SymmetryBoundaryConditionBlock(direction, side)
wall_const =  ["Minf","Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)

## Isothermal wall condition
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict = local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBoundaryConditionBlock(direction, 0, wall_eqns, local_dict)

# Top dirichlet shock condition
direction = 1; side = 1
pre_shock_pressure = 1.0/(1.4*4.0)
post_shock_pressure = pre_shock_pressure * 1.18646663
shock_loc = parse_expr("Eq(GridVariable(shock_loc), 40.0)", local_dict=local_dict)

d = "Eq(GridVariable(d), Piecewise((1.0 ,GridVariable(x0)<shock_loc),(1.129734572, GridVariable(x0)>=shock_loc),(1.0,True)))"
u0 = "Eq(GridVariable(u0), Piecewise((1.0 ,GridVariable(x0)<shock_loc),(0.9667023804242755, GridVariable(x0)>=shock_loc),(1.0,True)))"
u1 = "Eq(GridVariable(u1), Piecewise((0.0 ,GridVariable(x0)<shock_loc),(-0.052106102140246795, GridVariable(x0)>=shock_loc),(1.0,True)))"
p = "Eq(GridVariable(p), Piecewise( (%f,GridVariable(x0)<shock_loc),(%f, GridVariable(x0)>=shock_loc),(1.0,True)))" % (pre_shock_pressure, post_shock_pressure)

d = parse_expr(d, local_dict=local_dict); u0 = parse_expr(u0, local_dict=local_dict); u1 = parse_expr(u1, local_dict=local_dict)
p = parse_expr(p, local_dict=local_dict)

upper_eqns = [gridx0, gridx1, shock_loc, d, u0, u1, p, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBoundaryConditionBlock(direction, side, upper_eqns)

block.set_block_boundaries(boundaries)
# metric = copy.deepcopy(metriceq)
sim_eq = copy.deepcopy(simulation_eq)
CR = copy.deepcopy(constituent)

# block.set_equations([sim_eq, CR, metric])
block.set_equations([sim_eq, CR, initial])

block.set_discretisation_schemes(schemes)

block.discretise()

### Grid generation starts here 
# BL_profile = Boundary_layer_profile()
# [np0, np1], [L0, L1]

# Pass in block and required solution vector arrays to the grid gen
arrays = simulation_eq.time_advance_arrays

pprint(block.block_datasets)

init_grid = Initialise_Solution_On_Grid([400, 800], [350, 115], arrays, block)
# Direction 0, side 0, uniform, stretch factor 0 (no stretching)
# init_grid.apply_asymmetric_stretching(0, 0, 0)
# # Direction 1, side 0, stretch factor 5
# init_grid.apply_asymmetric_stretching(1, 0, 5.0)

init_grid.read_sbli_soln()
# init_grid.create_2D_meshgrid()
init_grid.add_coordinate_halos()


# init_grid.output_hdf5()


# Create boundary layer solution to interpolate onto the grid # 1.67619491
# BL = Katzer_setup(2.0, 0.72, 1.4, -1.0, 950.0, 2.0)
# BL.stretched_points = init_grid.stretched_points
# rho_arr, rhou_arr, rhov_arr, rhoE_arr = BL.setup()
# names = ['rho', 'rhou0', 'rhou1', 'rhoE']
# variables = [rho_arr, rhou_arr, rhov_arr, rhoE_arr]
# init_grid.set_solutions(names, variables)
init_grid.add_solution_halos()

init_grid.output_hdf5()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
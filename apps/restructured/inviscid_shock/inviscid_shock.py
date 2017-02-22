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


ndim = 2
weno_order = 5
weno = True
Euler_eq = EulerEquations(ndim, weno)
ev_dict, LEV_dict, REV_dict = Euler_eq.generate_eig_system()
Avg = SimpleAverage([0, 1])
LLF = LLFCharacteristic(ev_dict, LEV_dict, REV_dict, weno_order, ndim, Avg)

sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
a = "Conservative(rhou_j,x_j,%s)"%sc1
mass = "Eq(Der(rho,t), - %s)"%(a)
a = "Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s)"%sc1
momentum = "Eq(Der(rhou_i,t) , -%s  )"%(a)
a = "Conservative((p+rhoE)*u_j,x_j, %s)"%sc1
energy = "Eq(Der(rhoE,t), - %s  )"%(a)
# Substitutions
substitutions = []

# Define all the constants in the equations
constants = ["gama", "Minf"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"

simulation_eq = SimulationEquations()
eq = Equation()
eqns = eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(speed_of_sound, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)


schemes = {}
schemes[LLF.name] = LLF
rk = RungeKutta(3)
schemes[rk.name] = rk

block= SimulationBlock(ndim, block_number = 0)
block.sbli_rhs_discretisation = True

local_dict = {"block" : block, "GridVariable" : GridVariable, "DataObject" : DataObject}
for con in constants:
    local_dict[con] = ConstantObject(con)
# Initial conditions
x0 = parse_expr("Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
x1 = parse_expr("Eq(GridVariable(x1), block.deltas[1]*block.grid_indexes[1])", local_dict=local_dict)


rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)",local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5* d *(u0**2+u1**2))", local_dict=local_dict)



d_in = parse_expr("Eq(GridVariable(d), 1.0)", local_dict=local_dict)
u0_in = parse_expr("Eq(GridVariable(u0), 1.0)", local_dict=local_dict)
u1_in = parse_expr("Eq(GridVariable(u1), 0.0)", local_dict=local_dict)
p_in = parse_expr("Eq(GridVariable(p), 1.0/(gama*Minf**2.0))", local_dict=local_dict) # p = 1/(gama*Minf^2)

initial = GridBasedInitialisation()
initial.add_equations([d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE])

boundaries = [[0,0] for t in range(ndim)]
# Left Dirchlet at x= 0, inlet conditions
### Inflow conditions
direction = 0
side = 0
inlet_eq = [d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBoundaryConditionBlock(direction, side, inlet_eq)

# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = LinearExtrapolateBoundaryConditionBlock(direction, side)

# Bottom inviscid wall
direction = 1; side = 0 
boundaries[direction][side] = SymmetryBoundaryConditionBlock(direction, side)

# Top dirichlet shock condition
direction = 1; side = 1
pre_shock_pressure = 1.0/(1.4*4.0)
post_shock_pressure = pre_shock_pressure * 1.18646663
shock_loc = parse_expr("Eq(GridVariable(shock_loc), 40.0)", local_dict=local_dict)

d = "Eq(GridVariable(d), Piecewise((1.0 ,x0<shock_loc),(1.129734572, x0>=shock_loc),(1.0,True)))"
u0 = "Eq(GridVariable(u0), Piecewise((1.0 ,x0<shock_loc),(0.9667023804242755, x0>=shock_loc),(1.0,True)))"
u1 = "Eq(GridVariable(u1), Piecewise((0.0 ,x0<shock_loc),(-0.052106102140246795, x0>=shock_loc),(1.0,True)))"
p = "Eq(GridVariable(p), Piecewise( (%f,x0<shock_loc),(%f, x0>=shock_loc),(1.0,True)))" % (pre_shock_pressure, post_shock_pressure)

d = parse_expr(d, local_dict=local_dict); u0 = parse_expr(u0, local_dict=local_dict); u1 = parse_expr(u1, local_dict=local_dict)
p = parse_expr(p, local_dict=local_dict)

upper_eqns = [x0, shock_loc, d, u0, u1, p, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBoundaryConditionBlock(direction, side, upper_eqns)

block.set_block_boundaries(boundaries)
block.set_equations([copy.deepcopy(constituent),copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
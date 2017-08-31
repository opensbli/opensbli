#!/usr/bin/env python
import sys
import os
from math import ceil

# Import local utility functions
# import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.initialisation import *


BUILD_DIR = os.getcwd()
# Problem dimension
ndim = 3

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho, t), - Skew(rho*u_j, x_j))"
momentum = "Eq(Der(rhou_i, t), - Skew(rhou_i*u_j, x_j) - Der(p, x_i) + Der(tau_i_j, x_j) - KD(_i,_j)*c_j )"
energy = "Eq(Der(rhoE, t), - Skew(rhoE*u_j, x_j) - Conservative(p*u_j, x_j) - Dot(c_j, u_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j) )"
equations = [mass, momentum, energy]

stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

constants = ["Re", "Pr", "gama", "Minf", "mu", "c_j"]
coordinate_symbol = "x"

# metriceq = MetricsEquation()
# metriceq.genreate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 3)

velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(2/3))"
formulas = [velocity, pressure, temperature]

simulation_eq = SimulationEquations()
eq = Equation()
eqns = eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# simulation_eq.apply_metrics(metriceq)

# latex = LatexWriter()
# latex.open('./equation_transformations.tex')
# metadata = {"title": "Transformations of the equations in OpenSBLI framework", "author": "Satya P Jammy", "institution": "University of Southampton"}
# latex.write_header(metadata)
# for no, eq in enumerate(flatten(simulation_eq.equations)):
#    latex.write_expression(eq)
# latex.write_footer()
# latex.close()

constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

block = SimulationBlock(ndim, block_number=0)

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

x0 = "Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])"
x1 = "Eq(DataObject(x1), block.deltas[1]*block.grid_indexes[1])"
x2 = "Eq(DataObject(x2), block.deltas[2]*block.grid_indexes[2])"

# Laminar initial condition
# xl = 2.0*pi
# yl = 2.0
# zl = pi
x0l = "Eq(GridVariable(x0l), 2.0*pi)"
x1l = "Eq(GridVariable(x1l), 2.0)"
x2l = "Eq(GridVariable(x2l), pi)"

# u0 = "Eq(GridVariable(u0), 45*(1-(DataObject(x1)-1.0)**2))"
u0 = "Eq(GridVariable(u0), 0)"
u1 = "Eq(GridVariable(u1), 0)"
u2 = "Eq(GridVariable(u2), 0)"
p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf))"
r = "Eq(GridVariable(r), 1.0/(1.0+0.01944*(1-(DataObject(x1)-1)**4)))"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhou2 = "Eq(DataObject(rhou2), r*u2)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r*(u0**2 + u1**2 + u2**2))"
eqns = [x0, x1, x2, u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE]

initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
pprint(initial_equations)
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)


schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

boundaries = []
# Periodic boundaries in x direction
direction = 0
boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]

# Isothermal wall in y direction
direction = 1
rhoEd = "Eq(DataObject(rhoE), DataObject(rho)/((gama-1)*gama*Minf*Minf))"
rhoEd = parse_expr(rhoEd, local_dict=local_dict)
upper_wall_eq = [rhoEd]
lower_wall_eq = [rhoEd]
boundaries += [IsothermalWallBoundaryConditionBlock(direction, 0, upper_wall_eq, local_dict)]
boundaries += [IsothermalWallBoundaryConditionBlock(direction, 1, lower_wall_eq, local_dict)]

# Periodic boundaries in z direction
direction = 2
boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]

block.set_block_boundaries(boundaries)
# metric = copy.deepcopy(metriceq)
# block.set_equations([copy.deepcopy(constituent),copy.deepcopy(simulation_eq), metric, initial])
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

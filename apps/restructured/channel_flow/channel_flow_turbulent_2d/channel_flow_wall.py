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
ndim = 2

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

# Turbulent initial condition
# xl = 2.0*pi
# yl = 2.0
x0l = "Eq(GridVariable(x0l), 2.0*pi)"
x1l = "Eq(GridVariable(x1l), 2.0)"

sx0 = "Eq(GridVariable(sx0), sin(4.0*pi*DataObject(x0)/x0l) )"
sx1 = "Eq(GridVariable(sx1), sin(pi*(DataObject(x1)-1.0)) )"
cx0 = "Eq(GridVariable(cx0), cos(4.0*pi*DataObject(x0)/x0l) )"
cx1 = "Eq(GridVariable(cx1), 1.0+cos(pi*(DataObject(x1)-1.0)) )"
x1wall = "Eq(GridVariable(x1wall), Abs(1.0 - Abs(DataObject(x1)-1.0)) )"


vonkar = "Eq(GridVariable(vonkar), 2.5)"
b = "Eq(GridVariable(b), 5.5)"
visc = "Eq(GridVariable(visc), 1.0/180.0)"
amp = "Eq(GridVariable(amp), 0.1*(vonkar*log(1.0/visc)+b))"

ubar = "Eq(GridVariable(ubar), Piecewise((x1wall/visc, x1wall/visc < 10.0), (vonkar*log(x1wall/visc)+b, True)) )"

u0 = "Eq(GridVariable(u0), ubar+amp*x0l/2.0*cx0*sx1)"
u1 = "Eq(GridVariable(u1), -amp*sx0*sx1)"
p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r*(u0**2 + u1**2))"
eqns = [x0, x1, x0l, x1l, sx0, sx1, cx0, cx1, x1wall, vonkar, b, visc, amp, ubar, u0, u1, p, r, rho, rhou0, rhou1, rhoE]

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

block.set_block_boundaries(boundaries)
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
import copy
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.core.latex import *


BUILD_DIR = os.getcwd()

#base.LOG.info("Generating code for the 3D Taylor-Green Vortex simulation...")
#start_total = time.time()

# Problem dimension
ndim = 2
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.genreate_transformations(ndim, coordinate_symbol, [(False, False), (True, False)], 2)

mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j) )"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"

stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

substitutions = [stress_tensor, heat_flux]

constants = ["Re", "Pr","gama", "Minf", "mu"]

# Create a simulation equation class
simulation_eq = SimulationEquations()
# Add mass, momentum and Energy
simulation_eq.add_equations(Equation().expand(mass, ndim, coordinate_symbol, substitutions, constants))
simulation_eq.add_equations(Equation().expand(momentum, ndim, coordinate_symbol, substitutions, constants))
simulation_eq.add_equations(Equation().expand(energy, ndim, coordinate_symbol, substitutions, constants))
latex = LatexWriter()
original = flatten(simulation_eq.equations)
latex.open('./equation_transformations.tex')
metadata = {"title": "Transformations of the equations in OpenSBLI framework", "author": "Satya P Jammy", "institution": "University of Southampton"}
latex.write_header(metadata)
#simulation_eq.apply_metrics(metriceq)
for no, eq in enumerate(flatten(simulation_eq.equations)):
    latex.write_expression(original[no])
    latex.write_string(srepr(original[no]))
    latex.write_expression(eq)
latex.write_footer()
latex.close()
#exit()

velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"

# Create a constituent relations class
constituent = ConstituentRelations()
# Add the equations
constituent.add_equations(Equation().expand(velocity, ndim, coordinate_symbol, substitutions, constants))
constituent.add_equations(Equation().expand(pressure, ndim, coordinate_symbol, substitutions, constants))
constituent.add_equations(Equation().expand(temperature, ndim, coordinate_symbol, substitutions, constants))

# Create a simulation
block= SimulationBlock(ndim, block_number = 0)

boundaries = []
# Create boundaries, one for each side per dimension
for direction in range(ndim):
    boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
    boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]
# Set the boundary conditions in the x1 direction to adiabatic wall and carpenter by default
boundaries[2] = AdiabaticWall(1,0)
boundaries[3] = AdiabaticWall(1,1)
# Set the boundary conditions of the block
block.set_block_boundaries(boundaries)

# Create a Grid based initialisation object
initial = GridBasedInitialisation()
# Equations for the initial conditions
local_dict = {"block" : block, "GridVariable" : GridVariable, "DataObject" : DataObject}

x0 = "Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])"
x1 = "Eq(GridVariable(x1), block.deltas[1]*block.grid_indexes[1])"
x2 = "Eq(GridVariable(x2), block.deltas[2]*block.grid_indexes[2])"

u0 = "Eq(GridVariable(u0),sin(x0)*cos(x1)*cos(x2))"
u1 = "Eq(GridVariable(u1),-cos(x0)*sin(x1)*cos(x2))"
u2 = "Eq(GridVariable(u2), 0.0)"
p = "Eq(GridVariable(p), 1.0/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*x0)+cos(2.0*x1))*(2.0 + cos(2.0*x2)))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhou2 = "Eq(DataObject(rhou2), r*u2)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2 + u2**2))"

eqns = [x0, x1, x2, u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2, rhoE]
#initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
# Add the equations to the initialisation object
#initial.add_equations(initial_equations)

# set the equations to be solved on the block;
# We deepcopy the equations so that the original equations can be reused across the blocks
metric = copy.deepcopy(metriceq); simulation = copy.deepcopy(simulation_eq); CR = copy.deepcopy(constituent)
block.set_equations([simulation , CR, metric])
# Set the discretisation schemes for the block (A dictionary)
schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

block.set_discretisation_schemes(schemes)
# Discretise the equations on block
block.discretise()

# Algorithm for the block
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
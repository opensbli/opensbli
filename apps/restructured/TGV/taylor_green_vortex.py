#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *
from opensbli.initialisation import *


BUILD_DIR = os.getcwd()
# Problem dimension
ndim = 3
sc1 = "**{\'scheme\':\'Weno\'}"
sc2 = "**{\'scheme\':\'Central\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j,%s))" % sc2
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j, %s) - Der(p,x_i, %s)  + Der(tau_i_j,x_j,%s) )" % (sc2, sc2, sc2)


energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j,%s) - Conservative(p*u_j,x_j,%s) + Der(q_j,x_j,%s) + Der(u_i*tau_i_j ,x_j, %s) )" % (sc2, sc2, sc2, sc2)
ke = "Eq(ke, rho*(1/2)*Dot(u_j,u_j))"
enstrophy = "Eq(enstrophy, (1/2)*rho*(LC(_i,_j,_k)*Der(u_k,x_j,%s))**2)" % sc2
rhomean = "Eq(rhomean, rho)"
equations = [mass, momentum, energy]
diagnostics = [ke, enstrophy, rhomean]



stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j, %s)+ Der(u_j,x_i,%s)- (2/3)* KD(_i,_j)* Der(u_k,x_k,%s)))" % (sc2, sc2, sc2)
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j,%s))" % sc2

substitutions = [stress_tensor, heat_flux]

constants = ["Re", "Pr","gama", "Minf", "mu"]
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

constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)


# latex = LatexWriter()
# latex.open('./equations.tex')
# metadata = {"title": "Einstein Expansion of new Taylor green vortex equations", "author": "David", "institution": ""}
# latex.write_header(metadata)
# latex.write_string("Simulation equations\n")
# for index, eq in enumerate(flatten(simulation_eq.equations)):
#     if isinstance(eq, Equality):
#         latex.write_expression(eq)

# latex.write_string("Constituent relations\n")
# for index, eq in enumerate(flatten(constituent.equations)):
#     if isinstance(eq, Equality):
#         latex.write_expression(eq)
# latex.write_footer()
# latex.close()


block = SimulationBlock(ndim, block_number = 0)

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

initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
pprint(initial_equations)
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)


schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

block.sbli_rhs_discretisation = True

boundaries = []
# Create boundaries, one for each side per dimension
for direction in range(ndim):
	boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
	boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]

block.set_block_boundaries(boundaries)
block.set_equations([copy.deepcopy(constituent),copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)


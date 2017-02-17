#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *


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


latex = LatexWriter()
latex.open('./equations.tex')
metadata = {"title": "Einstein Expansion of new Taylor green vortex equations", "author": "David", "institution": ""}
latex.write_header(metadata)
latex.write_string("Simulation equations\n")
for index, eq in enumerate(flatten(simulation_eq.equations)):
    if isinstance(eq, Equality):
        latex.write_expression(eq)

latex.write_string("Constituent relations\n")
for index, eq in enumerate(flatten(constituent.equations)):
    if isinstance(eq, Equality):
        latex.write_expression(eq)
latex.write_footer()
latex.close()


block= SimulationBlock(ndim, block_number = 0)
# create Initial conditions
from opensbli.initialisation import *
initial = GridBasedInitialisation()
# As of now this is a work around later the initial conditions will be written as strings
xi = [GridVariable('x%d'%(i)) for i in range(block.ndim)] # 0 to ndim-1
le = [block.grid_indexes[i]*block.deltas[i] for i in range(block.ndim)]
xi += [GridVariable('u'), GridVariable('v'), GridVariable('w')] # Add the grid symbols index ndim to 2.*ndim-1
le += [sin(xi[0])*cos(xi[1])*cos(xi[2]), -cos(xi[0])*sin(xi[1])*cos(xi[2]), 0]
consts = [ConstantObject(r) for r in constants]
print consts
xi += [GridVariable('p')] # index in array 2*ndim
le += [1.0/(consts[2]*consts[3]**2) + (1.0/16.0)* ((cos(2.0*xi[0]) + cos(2.0*xi[1]))*(2.0 + cos(2.0*xi[2])))] # check this for 3d
xi += [GridVariable('r')] # index in array 2*ndim+1
le += [consts[2]*consts[3]**2*xi[2*ndim]]

arrs = flatten(simulation_eq.time_advance_arrays)

xi += [arrs[0]]; le += [xi[2*ndim+1]] # Density
xi += [arrs[1]]; le += [xi[2*ndim+1] * xi[ndim]] # uvelocity
xi += [arrs[2]]; le += [xi[2*ndim+1] * xi[ndim+1]] # vvelocity
xi += [arrs[3]]; le += [xi[2*ndim+1] * xi[ndim+2]] # wvelocity
xi += [arrs[4]]; le += [xi[2*ndim]/(consts[2]-1) + (1.0/2.0)*xi[2*ndim+1] *(xi[ndim]**2 + xi[ndim+1]**2 + xi[ndim+2]**2)  ] # Energy
eq1 = [Eq(x,y) for x,y in zip(xi, le)]
pprint(eq1)

initial.add_equations(eq1)

schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

block.sbli_rhs_discretisation = True
boundaries = [PeriodicBoundaryConditionBlock()]*2*ndim
block.set_block_boundaries(boundaries)
block.set_equations([copy.deepcopy(constituent),copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)


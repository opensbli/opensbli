#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *


BUILD_DIR = os.getcwd()

#base.LOG.info("Generating code for the 3D Taylor-Green Vortex simulation...")
#start_total = time.time()

# Problem dimension
ndim = 2
sc1 = "**{\'scheme\':\'Weno\'}"
sc2 = "**{\'scheme\':\'Central\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
a = "Conservative(rhou_j,x_j,%s)"%sc1
b = "Conservative(rhou_j,x_j,%s)"%sc2
mass = "Eq(Der(rho,t), - %s -%s)"%(a,b)
a = "Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s)"%sc1
b = "Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s)"%sc2
momentum = "Eq(Der(rhou_i,t) , -%s -%s  + Der(tau_i_j,x_j) )"%(a,b)
a = "Conservative((p+rhoE)*u_j,x_j, %s)"%sc1
b = "Conservative((p+rhoE)*u_j,x_j, %s)"%(sc2)
energy = "Eq(Der(rhoE,t), - %s - %s + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"%(a,b)
ke = "Eq(ke, rho*(1/2)*Dot(u_j,u_j))"
enstrophy = "Eq(enstrophy, (1/2)*rho*(LC(_i,_j,_k)*Der(u_k,x_j))**2)"
rhomean = "Eq(rhomean, rho)"
equations = [ momentum, energy]
diagnostics = [ke, enstrophy, rhomean]

# Substitutions
stress_tensor = "Eq(tau_i_j, ((mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k))))"
heat_flux = "Eq(q_j, (1/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr","gama", "Minf", "mu"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(2/3))"
#formulas = []
#metrics = [(True, True), (True, True), (False, False)]
#generate_metrics = MetricsEquation(3, coordinate_symbol, metrics, 2)
#exit()
simulation_eq = SimulationEquations()
eq = Equation()
eqns = eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)

simulation_eq.add_equations(eqns)

eqns = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
#exit()
#pprint(eqns)
simulation_eq.add_equations(eqns)
eqns = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
#exit()
#OpenSBLIExpression(simulation_eq.equations[0].rhs)
# Parse the constituent relations
constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
#pprint(constituent.equations)
#pprint(simulation_eq.equations)
# Discretise the constituent relations if it contains any derivatives
#simulation_eq.add_constituent_relations(constituent)
latex = LatexWriter()
latex.open('./equations.tex')
metadata = {"title": "Einstein Expansion of equations", "author": "Jammy", "institution": ""}
latex.write_header(metadata)
for index, eq in enumerate(flatten(simulation_eq.equations)):
    if isinstance(eq, Equality):
        latex.write_expression(eq)
latex.write_footer()
latex.close()
#exit()

block= SimulationBlock(ndim, block_number = 0)
# create Initial conditions
from opensbli.initialisation import *
initial = GridBasedInitialisation()
# As of now this is a work around later the initial conditions will be written as strings
xi = [GridVariable('x%d'%(i)) for i in range(block.ndim)] # 0 to ndim-1
le = [block.grid_indexes[i]*block.deltas[i] for i in range(block.ndim)]
xi += [GridVariable('u'), GridVariable('v')] # Add the grid symbols index ndim to 2.*ndim-1
le += [sin(xi[0])*cos(xi[1]), -cos(xi[0])*sin(xi[1])]
consts = [ConstantObject(r) for r in constants]
print consts
xi += [GridVariable('p')] # index in array 2*ndim
le += [1.0/(consts[2]*consts[3]**2) + (1.0/16.0)* (cos(2.0*xi[0]))] # check this for 3d
xi += [GridVariable('r')] # index in array 2*ndim+1
le += [consts[2]*consts[3]**2*xi[2*ndim]]

arrs = flatten(simulation_eq.time_advance_arrays)
xi += [arrs[0]]; le += [xi[2*ndim+1]] # Density
xi += [arrs[1]]; le += [xi[2*ndim+1] * xi[ndim]] # uvelocity
xi += [arrs[2]]; le += [xi[2*ndim+1] * xi[ndim+1]] # vvelocity
xi += [arrs[3]]; le += [xi[2*ndim]/(consts[2]-1) + (1.0/2.0)*xi[2*ndim+1] *(xi[ndim]**2 + xi[ndim+1]**2) ] # Energy
eq1 = [Eq(x,y) for x,y in zip(xi, le)]
# Equations for the constituent relations
#arrs = flatten(simulation_eq.time_advance_arrays)[1:]
#eq1 += [Eq(e, 0) for e in flatten(simulation_eq.time_advance_arrays)]
initial.add_equations(eq1)
#exit()
print isinstance(initial, NonSimulationEquations)
#copy.deepcopy(initial)
#exit()
## Create eigensystem
# Euler = EulerEquations(ndim)
# ev_dict, LEV_dict, REV_dict = Euler.generate_eig_system()
# weno_order = 5
# GLF = GLFCharacteristic(ev_dict, LEV_dict, REV_dict, weno_order)
# SF = ScalarLocalLFScheme(weno_order)
# # exit()

schemes = {}
cent = Central(4)
schemes[cent.name] = cent
# weno = Weno(weno_order)
# schemes[weno.name] = weno
rk = RungeKutta(3)
schemes[rk.name] = rk


weno_order = 3
weno = True
Euler_eq = EulerEquations(ndim, weno)
ev_dict, LEV_dict, REV_dict = Euler_eq.generate_eig_system() # Taking more time on my personal computer SPJ
#
Avg = SimpleAverage([0, 1])

LLF = LLFCharacteristic(ev_dict, LEV_dict, REV_dict, weno_order, ndim, Avg)

#
# Euler_eq.eq_to_vector_form(simulation_eq.equations)
# LLF.vector_notation = Euler_eq.vector_notation
# LLF.required_formulas = constituent.equations
schemes[LLF.name] = LLF

# TESTING generate weno scheme
block.sbli_rhs_discretisation = True
boundaries = [PeriodicBoundaryConditionBlock()]*2*ndim
block.set_block_boundaries(boundaries)
block.set_equations([copy.deepcopy(constituent),copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()
#exit()
#print 'block1 \n\n'
#pprint(simulation_eq.equations[0])
## This creates the traditional algorithm
#block1 = SimulationBlock(ndim, block_number = 1)
#block1.sbli_rhs_discretisation = True
#boundaries = [PeriodicBoundaryConditionBlock()]*2*ndim
#block1.set_block_boundaries(boundaries)
#block1.set_equations([copy.deepcopy(constituent),copy.deepcopy(simulation_eq)])
#block1.set_discretisation_schemes(schemes)
#block1.discretise()
#from .block import DataSetsToDeclare
print DataSetsToDeclare.datasetbases
#print block1.block_datasets
#print block.block_datasets
#for key in block.block_datasets:
    #print key, bl  ock.block_datasets[key], block.block_datasets[key].__dict__
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
# Create a descritisation class where the equations are deepcopied using deepcopy
#ex = copy.deepcopy(simulation_eq)
# block= SimulationBlock(ndim, block_number = 0)

print ConstantsToDeclare.constants
"""
Here I can do the algorithm on my own
[MainPrg, GridRead, Initialisation, Metriceq, Tloop, Start[Simulation_eq], RKloopst,
BC's, Simulatio_eq.spatial_kernels, Time_kernles], rkend, Diagnostics, output, tloopend,
output, diagnostics, programend]

"""
latex = LatexWriter()
latex.open('./equations.tex')
metadata = {"title": "Einstein Expansion of equations", "author": "Jammy", "institution": ""}
latex.write_header(metadata)
for index, eq in enumerate(flatten(simulation_eq.equations)):
    if isinstance(eq, Equality):
        latex.write_expression(eq)
latex.write_footer()
latex.close()


#################




#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *


BUILD_DIR = os.getcwd()

#base.LOG.info("Generating code for the 3D Taylor-Green Vortex simulation...")
#start_total = time.time()

# Problem dimension
ndim = 3

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,**{\'scheme\':\'Central\'}))"
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j)  + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
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
simulation_eq.add_equations(eqns)
eqns = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
#pprint(simulation_eq.equations)
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
#exit()
schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk
#print(ndim, type(ndim))
import copy
# Create a descritisation class where the equations are deepcopied using deepcopy
#ex = copy.deepcopy(simulation_eq)
block= SimulationBlock(ndim, block_number = 0)
block.sbli_rhs_discretisation = True
bounaries = [PeriodicBoundaryConditionBlock()]*2*ndim#, PeriodicBoundaryConditionBlock(), PeriodicBoundaryConditionBlock()]
block.set_block_boundaries(bounaries)
block.set_equations([constituent,simulation_eq])
block.set_discretisation_schemes(schemes)
block.discretise()
TraditionalAlgorithm(block)
latex = LatexWriter()
latex.open('./equations.tex')
metadata = {"title": "Einstein Expansion of equations", "author": "Jammy", "institution": ""}
latex.write_header(metadata)
for index, eq in enumerate(flatten(simulation_eq.equations)):
    if isinstance(eq, Equality):
        latex.write_expression(eq)
latex.write_footer()
latex.close()

v = "Eq(q22 ,wy(18)*((wy(11)*dxidx -wy(12)*detadx) +(-wy(8)*dxidy +wy(9)*detady)))"
v = parse_expr(v)
pprint(v)
print print_latex(v)

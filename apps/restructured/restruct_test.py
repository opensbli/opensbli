#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
import opensbli
from opensbli.problem import *
from opensbli.block import *

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 3D Taylor-Green Vortex simulation...")
start_total = time.time()

# Problem dimension
ndim = 3

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rho*u_j,x_j,**{\'scheme\':\'Central\'}) +  Conservative(rho*u_j,x_j,**{\'scheme\':\'Weno\'}))"
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j,x_j) - Der(p,x_i) + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t), - Conservative(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
ke = "Eq(ke, rho*(1/2)*Dot(u_j,u_j))"
enstrophy = "Eq(enstrophy, (1/2)*rho*(LC(_i,_j,_k)*Der(u_k,x_j))**2)"
rhomean = "Eq(rhomean, rho)"
equations = [ momentum, energy]
diagnostics = [ke, enstrophy, rhomean]

# Substitutions
stress_tensor = "Eq(tau_i_j, (1/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr","gama", "Minf"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(u_j*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(2/3))"
formulas = []
consta = {}
#pprint(local_dict)
#parsed = parse_expr(momentum , local_dict, standard_transformations, evaluate=False)
#pprint(parsed)
#pprint(srepr(parsed))
# Create the TGV problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, formulas)
#expanded_equations = problem.get_expanded(problem.equations)
#expanded_formulas = problem.get_expanded(problem.formulas)
#pprint(flatten(expanded_equations)[1].atoms(ConstantObject))
#pprint(flatten(expanded_equations)[1].atoms(EinsteinTerm))
#pprint(srepr(expanded_equations[0]))
#spat, temproa = get_derivatives(flatten(expanded_equations))
#pprint((spat[0].atoms(DataObject)))
#pprint(
#for d in spat+temproa:
    #v = d.subderivatives()
    #if v:
        #pprint([d,'subevals', v])
    #else:
        #pprint ([d,'No sub evaluations', v])
exit()
# Create blocboundary conditions
from opensbli.bcs import *

bounaries = [PeriodicBoundaryCondition(), PeriodicBoundaryCondition(),Wall(), SymmetryBoundaryCondition(), PeriodicBoundaryCondition(), PeriodicBoundaryCondition()]
block0 = SimulationBlock(ndim)
block0.set_block_number(0)
block0.set_block_boundaries(bounaries)
block0.set_constituent_equations(expanded_formulas)

block0.instantiate_block()

# Create the spatial scheme for local equations
from opensbli.central import *
#spatial_scheme = Central(4)
central_equations = [eq.rhs for eq in flatten(expanded_equations)]


from opensbli.timestepping import *
n = 1 # number of spatial schemes per block
group = SpatialDescritisationGroup(n)
group[0] = SpatialCentral(4)
group[0].set_equations(central_equations, block0)

n = 1 # number of temporal schemes per block this would be generally 1

temporalgroup = RungeKutta(3)
temporalgroup.set_derivatives([eq.lhs for eq in flatten(expanded_equations)])
# number of diagnostics per block
diagnostics = DiagnosticsCentral(4, central_equations)

block0.descritise(group, temporalgroup)
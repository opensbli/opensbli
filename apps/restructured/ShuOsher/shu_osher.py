#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from opensbli.core.weno_opensbli import *
import copy

ndim = 1
sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
a = "Conservative(rhou_j,x_j,%s)" % sc1
mass = "Eq(Der(rho,t), - %s)" % (a)
a = "Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s)" % sc1
momentum = "Eq(Der(rhou_i,t) , -%s  )" % (a)
a = "Conservative((p+rhoE)*u_j,x_j, %s)" % sc1
energy = "Eq(Der(rhoE,t), - %s  )" % (a)
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
eq = EinsteinEquation()
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


block = SimulationBlock(ndim, block_number=0)

# Initial conditions
initial = GridBasedInitialisation()
# x = "GridVariable(x0)"
x0 = "Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])"
p = "Eq(GridVariable(p), Piecewise((10.33, x0<1.0),(1.0, x0>=1.0),(0.0,True)))"
u0 = "Eq(GridVariable(u0), Piecewise((2.629369, x0<1.0),(0.0, x0>=1.0),(0.0,True)))"
d = "Eq(GridVariable(d), Piecewise((3.857143, x0<1.0),(1.0+0.2*sin(5*GridVariable(x0)), x0>=1.0),(0,True)))"
rho = "Eq(DataObject(rho), d)"
rhou0 = "Eq(DataObject(rhou0), d*u0)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1.0) + 0.5* d *(u0**2.0))"

eqns = [x0, u0, p, d, rho, rhou0, rhoE]

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)


# Shu Osher boundary condition values left side
arrays = flatten(simulation_eq.time_advance_arrays)
subs_dict = {Symbol('x0'): 0}
left_eqns = [eq.subs(subs_dict) for eq in initial_equations]
pprint(left_eqns)

subs_dict = {Symbol('x0'): 10.0}

right_eqns = [eq.subs(subs_dict) for eq in initial_equations]

pprint(right_eqns)

boundaries = []
# Create boundaries, one for each side per dimension
for direction in range(ndim):
    boundaries += [DirichletBoundaryConditionBlock(direction, 0, left_eqns)]
    boundaries += [DirichletBoundaryConditionBlock(direction, 1, right_eqns)]


schemes = {}
# Local LaxFredirich scheme for weno
weno_order = '5Z'
# Generate the Eigen system for the Euler equations
Euler_eq = EulerEquations(ndim)
ev_dict, LEV_dict, REV_dict = Euler_eq.generate_eig_system()
# Averaging procedure to be used for the eigen system evaluation
Avg = SimpleAverage([0, 1])
# LLF scheme
LLF = LLFWeno(ev_dict, LEV_dict, REV_dict, weno_order, ndim, Avg)
# Add to schemes
schemes[LLF.name] = LLF
rk = RungeKutta(3)
schemes[rk.name] = rk

block.set_block_boundaries(boundaries)
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

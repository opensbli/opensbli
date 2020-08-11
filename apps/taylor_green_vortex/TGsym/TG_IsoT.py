#!/usr/bin/env python

# This opensbli set-up uses an isothermal equation of state (i.e. no energy equation is solved), with the speed of sound set
# by the specified Mach number. The test case is the standard triply-symmetric Taylor-Green problem at Re=1600.
# NDS August 2020

# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Number of dimensions of the system to be solved
ndim = 3

# Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
# Specify the schemes
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j))"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"

# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"

#substitutions = [stress_tensor, heat_flux]
substitutions = [stress_tensor]

# Constants that are used
constants = ["Re", "gama", "Minf"]

# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p,rho/(gama*Minf*Minf))" # barotropic equation of state


# Instantiate EinsteinEquation class for expanding the Einstein indices in the equations
einstein_eq = EinsteinEquation()

# Expand the simulation equations, for this create a simulation equations class
simulation_eq = SimulationEquations()

# Expand mass and add the expanded equations to the simulation equations
eqns = einstein_eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand momentum add the expanded equations to the simulation equations
eqns = einstein_eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations()  # Instantiate constituent relations object

# Expand momentum add the expanded equations to the constituent relations
eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Write the expanded equations to a Latex file with a given name and titile
latex = LatexWriter()
latex.open('equations.tex', "Einstein Expansion of the simulation equations")
latex.write_string("Simulation equations\n")
for index, eq in enumerate(flatten(simulation_eq.equations)):
    latex.write_expression(eq)

latex.write_string("Constituent relations\n")
for index, eq in enumerate(flatten(constituent.equations)):
    latex.write_expression(eq)

latex.close()

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
x0 = "Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])"
x1 = "Eq(GridVariable(x1), block.deltas[1]*block.grid_indexes[1])"
x2 = "Eq(GridVariable(x2), block.deltas[2]*block.grid_indexes[2])"

u0 = "Eq(GridVariable(u0),sin(x0)*cos(x1)*cos(x2))"
u1 = "Eq(GridVariable(u1),-cos(x0)*sin(x1)*cos(x2))"
u2 = "Eq(GridVariable(u2), 0.0)"
p = "Eq(GridVariable(p), 1/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*x0)+cos(2.0*x1))*(2.0 + cos(2.0*x2)))"
r = "Eq(GridVariable(r), gama*p*Minf*Minf)"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhou2 = "Eq(DataObject(rhou2), r*u2)"

eqns = [x0, x1, x2, u0, u1, u2, p, r, rho, rhou0, rhou1, rhou2]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKutta(3)
schemes[rk.name] = rk

boundaries = []
# Create boundaries, one for each side per dimension, so in total 6 BC's for 3D'
for direction in range(ndim):
    boundaries += [SymmetryBC(direction, 0)]
    boundaries += [SymmetryBC(direction, 1)]

# set the boundaries for the block
block.set_block_boundaries(boundaries)

# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=100, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
# h5.add_arrays([DataObject('x0'), DataObject('x1'), DataObject('x2')])
block.setio(copy.deepcopy(h5))
# set the equations to be solved on the block
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)

constants = ['Re', 'gama', 'Minf', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0']
values = ['1600.0', '1.4', '0.5', '0.02', '1000', '49', '49', '49', 'M_PI/(block0np0-1)', 'M_PI/(block0np1-1)', 'M_PI/(block0np2-1)']
substitute_simulation_parameters(constants, values)

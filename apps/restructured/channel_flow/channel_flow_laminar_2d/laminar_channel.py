#!/usr/bin/env python

# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# Problem dimension
ndim = 2

# Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
# Specify the schemes
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j))"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j) - KD(_i,_j)*c_j)"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j)  - Dot(c_j, u_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"

# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

substitutions = [stress_tensor, heat_flux]

# Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "mu", "c_j"]

# symbol for the coordinate system in the equations 
coordinate_symbol = "x"

# Constituent relations used in the system
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"

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

# Expand energy equation add the expanded equations to the simulation equations
eqns = einstein_eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

# Expand the constituent relations and them to the constituent relations class
constituent = ConstituentRelations() # Instantiate constituent relations object

# Expand momentum add the expanded equations to the constituent relations
eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
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
# Define the variables used for creating boundary conditions and the initialisation
# dx and dy of the grid
dx, dy = block.deltas
# Indices for the grid location
i,j = block.grid_indexes
# Some constants used 
gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})
""" Conservative vector is the time advancement arrays of the simulation equations.
the order follows the order in which they are added to the simulation equations
class, i.e. arrays of density, momentum (components), energy in the present case 
"""
q_vector = flatten(simulation_eq.time_advance_arrays)

# STEP 1
# Set the boundary conditions on the block
boundaries = []
# For laminar channel flow case the boundaries are periodic in x and walls in y
# Periodic boundaries in x0 direction
direction = 0
boundaries += [PeriodicBC(direction, side=0)]
boundaries += [PeriodicBC(direction, side=1)]

# Isothermal wall in x1 direction
# Energy on the wall is set 
wall_energy = [Eq(q_vector[3], q_vector[0]/(gama * Minf**2.0 * (gama - S.One)))]

direction = 1
# Side 0 (bottom wall) boundary
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]

# Side 1 (top) boundary
upper_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 1, upper_wall_eq)]

# set the boundaries for the block
block.set_block_boundaries(boundaries)

# The equation classes used for the block, these are
# simulation equations, constituent relations, grid and initial conditions 

# Create the grid and intial conditions
# Arrays to store x and y coordinates, i.e (x0 and x1)
x,y = symbols('x0:%d' %ndim, **{'cls':DataObject})
grid_equations = []
# Equations for generating the grid, simple equispacing grid
grid_equations += [Eq(x, i*dx), Eq(y, j*dy)]

# Initialisation equations
initial_equations = []
# local varibales for temperature and pressure
temperature, pressure = symbols('T p', **{'cls':GridVariable})
# Equations for pressure and temperature
initial_equations += [Eq(pressure, S.One/(gama*Minf**2.0))]
initial_equations += [Eq(temperature, S.One + 0.01944*(S.One - (y - S.One)**4))]

# Initialise the conservative vector
initial_equations += [Eq(q_vector[0], S.One/temperature)]
initial_equations += [Eq(q_vector[1], S.Zero)]
initial_equations += [Eq(q_vector[2], S.Zero)]
initial_equations += [Eq(q_vector[3], pressure/(gama - S.One))]

# Instantiate a grid based initialisation classes
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)

# STEP 2
# Set the equation classes for the block (list)
block.set_equations([constituent, simulation_eq, initial])


# STEP 3
# Create the dictionary of schemes
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKutta(3)
schemes[rk.name] = rk
# Set the discretisation schemes to be used (a python dictionary)
block.set_discretisation_schemes(schemes)

# STEP 4 add io for the block
kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays + [x,y]
output_hdf5 = iohdf5(arrays=output_arrays, **kwargs)
block.setio([output_hdf5])

# STEP 6
# Perform the symbolic discretisation of the equations 
block.discretise()

# STEP 7
# create an algorithm from the numerical solution
alg = TraditionalAlgorithmRK(block)

# STEP 8 
# set the simulation data type: if not set "Double" is default
SimulationDataType.set_datatype(Double)

# STEP 9
# Write the OPSC compatible code for the numerical solution
OPSC(alg)

# STEP 10 
# Populate the values of the constants like Re, Pr etc and the number of points for the
# simulation etc. In the future reading thes from HDF5 would be provided

from opensbli.utilities.helperfunctions import substitute_simulation_parameters
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', "c0", "c1"]
values = ['90.0', '1.4', '0.01', '0.72', '0.0001', '3000000', '16', '64', '2.0*M_PI/block0np0', '2.0/(block0np1-1)', '-1', '0']
substitute_simulation_parameters(constants, values)
exit()

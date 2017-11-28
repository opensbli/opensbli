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
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"

# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

substitutions = [stress_tensor, heat_flux]

# Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "mu"]

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

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
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

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
pprint(initial_equations)
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

# Boundary conditions
boundaries = []
# Periodic boundaries in x0 direction
direction = 0
boundaries += [PeriodicBC(direction, 0)]
boundaries += [PeriodicBC(direction, 1)]

# Isothermal wall in x1 direction
direction = 1
rhoEd = "Eq(DataObject(rhoE), DataObject(rho)/((gama-1)*gama*Minf*Minf))"
rhoEd = parse_expr(rhoEd, local_dict=local_dict)
upper_wall_eq = [rhoEd]
lower_wall_eq = [rhoEd]
boundaries += [IsothermalWallBC(direction, 0, upper_wall_eq)]
boundaries += [IsothermalWallBC(direction, 1, lower_wall_eq)]

# set the boundaries for the block
block.set_block_boundaries(boundaries)
# set the equations to be solved on the block
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])
# set the discretisation schemes 
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

# create an algorithm from the discretised computations, this writes the file algorithm.tex in latex_output directory
alg = TraditionalAlgorithmRK(block)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)

constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', "c0", "c1"]
values = ['180.0', '1.4', '0.01', '0.72', '0.00001', '1000000', '1024', '1024', '2.0*M_PI/block0np0', '2.0/(block0np1-1)', '-1', '0']
substitute_simulation_parameters(constants, values)

#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import pi, sin, cos, Abs, sqrt

# Problem dimension
ndim = 2
# Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central
mass     = "Eq(Der(rho,t), -(rhob*Der(u_i,x_i) + rho*Der(ub_i,x_i) + ub_i*Der(rho,x_i) + u_i*Der(rhob,x_i)) )"
momentum = "Eq(rhob*Der(u_i,t), -(rhob*u_j*Der(ub_i,x_j) + rho*ub_j*Der(ub_i,x_j) + (1.0/(gama*Minf**2))*( rhob*Der(T,x_i) + rho*Der(Tb,x_i) + Tb*Der(rho,x_i) + T*Der(rhob,x_i) ) - (1.0/Re)*Der(tau_i_j,x_j)))"
energy   = "Eq(rhob*Der(T,t), -(rhob*u_j*Der(Tb,x_j) + rho*ub_j*Der(Tb,x_j) + (gama-1)*(rhob*T*Der(ub_j,x_j) + rho*Tb*Der(ub_j,x_j)) - ((Minf**2*gama*(gama-1))/Re)*(taub_i_j*Der(u_i,x_j) + tau_i_j*Der(ub_i,x_j)) + (gama/(Re*Pr))*(Der(mub*Der(T,x_j),x_j) + dm_j*Der(Tb,x_j) + dmubTb*T*Der(Der(Tb,x_k),x_k))))"

# Substitutions used in the equations
stress_tensorb = "Eq(taub_i_j, (Der(ub_i,x_j)+ Der(ub_j,x_i)- (2/3)* KD(_i,_j)* Der(ub_k,x_k)))"
stress_tensor  = "Eq(tau_i_j, mub*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)) + dmubTb*T*(Der(ub_i,x_j)+ Der(ub_j,x_i)- (2/3)* KD(_i,_j)* Der(ub_k,x_k)) )"

substitutions = [stress_tensorb, stress_tensor]

# Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]

# symbol for the coordinate system in the equations
coordinate_symbol = "x"

# Constituent relations used in the system
# Tb         = "Eq(Tb, 1)" # The base flow temperature is now evaluated in the initial condition
viscosityb = "Eq(mub, Tb**(1.5)*(1.0+SuthT/RefT)/(Tb+SuthT/RefT))"
dmubTb     = "Eq(dmubTb, 1.5*Tb**(0.5)*(1.0+SuthT/RefT)/(Tb+SuthT/RefT) - Tb**(1.5)*(1.0+SuthT/RefT)/(Tb+SuthT/RefT)**2)"
# d2mubTb    = "Eq(d2mubTb,  )" # I've copied the evaluation of d2mubTb directly into the dmi evaluations
# viscosity  = "Eq(mu,  dmubTb*T)" # I've replaced mu with dmubTb*T in the governing equations instead. Is the 1.0/Re correct in the equations?

## This was called dmuxi before
dmi      = "Eq(dm_i, dmubTb*Der(T,x_i) + (3.0/4.0)*((1.0+SuthT/RefT)/(Tb**0.5*(Tb+SuthT/RefT))) - 3.0*((Tb**0.5*(1.0+SuthT/RefT))/((Tb+SuthT/RefT)**2)) + 2.0*((Tb**1.5*(1.0+SuthT/RefT))/((Tb+SuthT/RefT)**3))*Der(Tb,x_i)*T)"

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
constituent = ConstituentRelations()  # Instantiate constituent relations object

# Expand viscosity (base flow) add the expanded equations to the constituent relations
# eqns = einstein_eq.expand(Tb, ndim, coordinate_symbol, substitutions, constants)
# constituent.add_equations(eqns)

eqns = einstein_eq.expand(viscosityb, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = einstein_eq.expand(dmubTb, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# eqns = einstein_eq.expand(d2mubTb, ndim, coordinate_symbol, substitutions, constants)
# constituent.add_equations(eqns)

# eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, substitutions, constants)
# constituent.add_equations(eqns)

eqns = einstein_eq.expand(dmi, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
simulation_eq.apply_metrics(metriceq)

for eqn in constituent.equations:
    pprint(eqn)

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
u0 = "Eq(DataObject(u0),0.0)"
u1 = "Eq(DataObject(u1),0.0)"
T  = "Eq(DataObject(T), 0.0)"
r  = "Eq(DataObject(rho), 0.0)"

ub0  = "Eq(DataObject(ub0),  DataObject(rhou0)/DataObject(rho))"
ub1  = "Eq(DataObject(ub1),  DataObject(rhou1)/DataObject(rho))"
Tb   = "Eq(DataObject(Tb),   (gama*Minf*Minf*(gama-1)*(DataObject(rhoE)-0.5*DataObject(rho)*(DataObject(ub0)**2 + DataObject(ub1)**2)))/DataObject(rho))"
rhob = "Eq(DataObject(rhob), DataObject(rho))"

# The list ordering determines the order these equations will be evaluated in
eqns = [ub0, ub1, Tb, rhob, u0, u1, T, r] # Set the base flow first to avoid read/write conflicts with rho

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
## Low storage optimisation for the central scheme
fns = 'u0 u1 ub0 ub1 T Tb'
# cent = StoreSome(4, fns)
cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension
q_vector = flatten(simulation_eq.time_advance_arrays)
boundaries = []
# Apply a periodic boundary over the shared mesh line
boundaries += [PeriodicBC(direction=0, side=0)]
boundaries += [PeriodicBC(direction=0, side=1)]
# Isothermal wall in x1 direction
boundaries += [PrimitiveIsothermalWallBC(direction=1, side=0, wall_temperature=1.0)]
# Far field boundary
dirichlet_equations = [Eq(GridVariable('temp'), DataObject('x0'))]
boundaries += [DirichletBC(direction=1, side=1, equations=dirichlet_equations)]
# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=1000000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1'),  DataObject('rho'), DataObject('u0'), DataObject('u1'), DataObject('T')])
kwargs = {'iotype': "Read"}
h5_read = iohdf5(filename='base_flow.h5', **kwargs)
h5_read.add_arrays([DataObject('x0'), DataObject('x1'),DataObject('rho'),DataObject('rhou0'),DataObject('rhou1'),DataObject('rhoE')])
block.setio([h5, h5_read])
# Set the equations to be solved on the block
block.set_equations([constituent, simulation_eq, initial, metriceq])
# set the discretisation schemes
block.set_discretisation_schemes(schemes)
# Discretise the equations on the block
block.discretise()

arrays = ['u1', 'u1', 'u1', 'u1', 'u1', 'u1', 'u1', 'u1']
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(178, 44), (178, 73), (178, 97), (178, 121), (178, 143), (178, 165), (178, 187), (178, 198) ]
SM = SimulationMonitor(arrays, indices, block, print_frequency=250, fp_precision=12, output_file='output.log')
# create an algorithm from the discretised computations
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)

# Simulation parameters
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'SuthT', 'RefT']
values = ['100.0', '1.4', '0.1', '0.71', '0.0001', '10000000', '357', '201', '120.0/(block0np0-1)', '120.0/(block0np1-1)', '110.4', '288.0']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho_B0')

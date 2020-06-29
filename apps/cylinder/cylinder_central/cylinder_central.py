#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import pi, sin, cos, Abs, sqrt

# Sponge boundary condition on the inflow
# def Sponge_Boundary_Kernel(conservative_variables, block):
#     sponge_length, Lc, damping_factor, r_coord = symbols("Ls Lc damping_factor r", **{'cls':GridVariable})
#     gama, Minf = symbols("gama Minf", **{'cls':ConstantObject})
#     Residual = symbols("Residual0:4", **{'cls':DataObject})
#     # Set the length of the sponge zone
#     equations = [OpenSBLIEq(sponge_length, 5.0)]
#     # Calculate the radial length
#     equations += [OpenSBLIEq(r_coord, sqrt(DataObject('x0')**2 + DataObject('x1')**2))]
#     equations += [OpenSBLIEq(Lc, r_coord - (60.0 - sponge_length))]
#     # Damping equation
#     equations += [OpenSBLIEq(damping_factor,  0.5*(1.0 + cos(pi * Lc/sponge_length)))]
#     # Prescribed inflow conditions (rho, rhou, rhov, rhoE)
#     freestream_values = [1.0, 1.0, 0.0, 1.0/(gama * Minf**2.0 *(gama - 1.0)) + 0.5]
#     # Calculate the correction if within the sponge zone
#     i = 0
#     for R, Q, Q_star in zip(Residual, conservative_variables, freestream_values):
#     	condition1, condition2 = ExprCondPair(damping_factor * (Q - Q_star), r_coord >= 55.0), ExprCondPair(0.0, True)
#         equations += [OpenSBLIEq(GridVariable('correction_%d' % i), Piecewise(condition1, condition2))]
#         i += 1
#     # Apply the correction
#     for i, R in enumerate(Residual):
#     	equations += [OpenSBLIEq(R, R - GridVariable('correction_%d' % i))]
#     # Create a kernel and set the equations
#     eqns = block.dataobjects_to_datasets_on_block(equations)
#     ker = Kernel(block, computation_name="Sponge kernel")
#     ker.kernelname = "sponge_boundary_condition"
#     ker.add_equation(eqns)
#     ranges = copy.deepcopy(block.ranges)
#     ker.ranges = ranges
#     # for eqn in ker.equations:
#     # 	pprint(eqn)
#     ker.update_block_datasets(block)
#     return ker

# Problem dimension
ndim = 2
# Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
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
constituent = ConstituentRelations()  # Instantiate constituent relations object
# Expand momentum add the expanded equations to the constituent relations
eqns = einstein_eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
# Expand pressure add the expanded equations to the constituent relations
eqns = einstein_eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
# Expand temperature add the expanded equations to the constituent relations
eqns = einstein_eq.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)
simulation_eq.apply_metrics(metriceq)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)

# Local dictionary for parsing the expressions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

# Initial conditions as strings
u0 = "Eq(GridVariable(u0),1.0)"
u1 = "Eq(GridVariable(u1), 0.0,)"
p = "Eq(GridVariable(p), 1/(gama*Minf*Minf))"
r = "Eq(GridVariable(r), gama*Minf*Minf*p)"

rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r *(u0**2+ u1**2))"

eqns = [u0, u1, p, r, rho, rhou0, rhou1, rhoE]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(initial_equations)

# Create a schemes dictionary to be used for discretisation
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
# Low storage optimisation for the central scheme
fns = 'u0 u1 T'
cent = StoreSome(4, fns)
# cent = Central(4)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk

# Create boundaries, one for each side per dimension
q_vector = flatten(simulation_eq.time_advance_arrays)
boundaries = []
direction = 0
# Apply a periodic boundary over the shared mesh line
boundaries += [PeriodicBC(direction, 0)]
boundaries += [PeriodicBC(direction, 1)]
# Isothermal wall in x1 direction
gama, Minf = symbols('gama Minf', **{'cls': ConstantObject})
# Energy on the wall is set
wall_energy = [Eq(q_vector[3], q_vector[0] / (gama * Minf**2.0 * (gama - S.One)))]
direction = 1
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]
# Far field boundary
direction, side = 1,1
boundaries += [DirichletBC(direction, side, initial_equations)]
# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=100000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')]) #,  DataObject('rho_filt'), DataObject('rhou0_filt'), DataObject('rhou1_filt'), DataObject('rhoE_filt')])
kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio([h5, h5_read])

# Add SFD filtering
# SFD = SFD(block, chifilt=0.1, omegafilt=1.0/0.75)
j = block.grid_indexes[1]
grid_condition = j >= 169
F = BinomialFilter(block, order=10, grid_condition=grid_condition)

# Set the equations to be solved on the block
block.set_equations([constituent, simulation_eq, initial, metriceq] + F.equation_classes)# + SFD.equation_classes)
# set the discretisation schemes
block.set_discretisation_schemes(schemes)

# Discretise the equations on the block
block.discretise()

# Add a kernel for the sponge boundary condition after block.discretise()
# Q = flatten(simulation_eq.time_advance_arrays)
# sponge_kernel = Sponge_Boundary_Kernel(Q, block)
# for no, eq_class in enumerate(block.list_of_equation_classes):
#     if isinstance(eq_class, SimulationEquations):
#         eq_class.Kernels += [sponge_kernel]

arrays = ['p', 'p', 'p', 'p', 'p', 'p', 'p']
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(178, 45), (178, 72), (178, 96), (178, 118), (178, 139), (178, 160), (178, 176)]
SM = SimulationMonitor(arrays, indices, block, print_frequency=250, fp_precision=12, output_file='output.log')
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# set the simulation data type, for more information on the datatypes see opensbli.core.datatypes
SimulationDataType.set_datatype(Double)

# Write the code for the algorithm
OPSC(alg)
# Change grid size here if desired
f = h5py.File('grid.h5', 'r')
x0, x1 = f['x0'].value, f['x1'].value
print(x1.shape)
npoints = [357, 179]
halos = [(-5, 5), (-5, 5)]
arrays, array_names = [x0, x1], ['x0', 'x1']
output_hdf5(arrays, array_names, halos, npoints, block)
# Simulation parameters
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0']
values = ['100.0', '1.4', '0.1', '0.71', '0.0001', '1000000', '357', '179', '120.0/(block0np0-1)', '120.0/(block0np1-1)']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho_B0')

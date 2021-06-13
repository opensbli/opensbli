#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from sympy import sin, log, cos, pi
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

# STEP 0 Create the equations required for the numerical solution
# Problem dimension
ndim = 3
stats = True

# Define the compresible Navier-Stokes equations in Einstein notation
# Feiereisen quadratic skew-symmetric formulation, no change in continuity
mass = "Eq(Der(rho, t), - Der(rhou_j, x_j))"

# Feiereisen quadratic skew-symmetric momentum
# TODO add the refernece paper
# we expand convective and viscous parts separately and add them later
# this demonstrates how OpenSBLI equations can be used to build equations

QSSFm = "(1/2) * (Conservative(rhou_i*u_j, x_j) + rhou_j* Der(u_i,x_j) + u_i * Der(rhou_j,x_j))"
momentum = "Eq(Der(rhou_i, t), - Der(p, x_i) + Der(tau_i_j, x_j) - KD(_i,_j)*c_j )"

QSSFe = "(1/2) * (Conservative(rhoE*u_j, x_j) + rhou_j*Conservative(rhoE/rho, x_j) + (rhoE/rho) * Der(rhou_j, x_j))"
energy = "Eq(Der(rhoE, t), - %s - Conservative(p*u_j, x_j) - Dot(c_j, u_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j) )" % (QSSFe)

stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"

substitutions = [stress_tensor, heat_flux]
# Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "c_j"]

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
expanded_Feiereisen = einstein_eq.expand(QSSFm, ndim, coordinate_symbol, substitutions, constants)
eqns = einstein_eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
# Substract the inviscid part ot the RHS
for no, value in enumerate(eqns):
    eqns[no] = Eq(eqns[no].lhs,  eqns[no].rhs - expanded_Feiereisen[no])
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

# Write the expanded equations to a Latex file with a given name and titile
latex = LatexWriter()
latex.open('equations.tex', "Einstein Expansion of the equations")
simulation_eq.write_latex(latex)
constituent.write_latex(latex)
latex.close()

if stats:
    # Create the statistics equations, this shows another way of writing the equations
    from stats import *
    q_vector = flatten(simulation_eq.time_advance_arrays)

    # Create equation classes for statistics, see stats.py in the current folder.
    # Later an automatic way of creating statistics equations will be provided
    stat_equation_classes = favre_averaged_stats(ndim, q_vector)
else:
    stat_equation_classes = []
# Create a simulation block
block = SimulationBlock(ndim, block_number=0)

# Define the variables used for creating boundary conditions and the initialisation
# dx and dy of the grid
dx, dy, dz = block.deltas
# Indices for the grid location
i, j, k = block.grid_indexes
# Some constants used
gama, Minf, Re = symbols('gama Minf Re', **{'cls': ConstantObject})
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
Twall = ConstantObject("Twall")
wall_energy = [Eq(q_vector[4], Twall*q_vector[0] / (gama * Minf**2.0 * (gama - S.One)))]

direction = 1
# Side 0 (bottom wall) boundary
lower_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 0, lower_wall_eq)]

# Side 1 (top) boundary
upper_wall_eq = wall_energy[:]
boundaries += [IsothermalWallBC(direction, 1, upper_wall_eq)]

# Periodic boundaries in x2 direction
direction = 2
boundaries += [PeriodicBC(direction, 0)]
boundaries += [PeriodicBC(direction, 1)]

# set the boundaries for the block
block.set_block_boundaries(boundaries)

# Create the grid and intial conditions
# Arrays to store x and y coordinates, i.e (x0 and x1)
x, y, z = symbols('x0:%d' % ndim, **{'cls': DataObject})
grid_equations = []
# Equations for generating the grid, simple equispacing grid, later change to stretched
grid_equations += [Eq(x, i * dx), Eq(y, -S.One + j * dy), Eq(z, k * dz)]

# Initial condition
initial_equations = []

lx, ly, lz = symbols('lx0:%s' % ndim, **{'cls': ConstantObject})
sx, sy, sz = symbols('sx0:%s' % ndim, **{'cls': GridVariable})
cx, cy, cz = symbols('cx0:%s' % ndim, **{'cls': GridVariable})
ubar, amp, vonkar, b = symbols('ubar amp vonkar b', **{'cls': GridVariable})

# Set the constants for initialisation
initial_equations += [Eq(b, 5.5), Eq(vonkar, 2.5)]

# Equation for umean
initial_equations += [Eq(ubar, Piecewise((y * Re, y * Re < 10.0),
                         (vonkar * log(y * Re) + b, True)))]
# Amplitude of the disturbances
initial_equations += [Eq(amp, 0.1 * (vonkar * log(Re) + b))]

# sin disturbances in each direction
initial_equations += [Eq(sx, sin(4.0 * pi * x / lx)), Eq(sy, sin(pi * y)),
                         Eq(sz, sin(2.0 * pi * z / lz))]

# cosine disturbances in each direction
initial_equations += [Eq(cx, cos(4.0 * pi * x / lx)), Eq(cy, S.One + \
                         cos(pi * y)), Eq(cz, cos(2.0 * pi * z / lz))]

d, u, v, w, p = symbols('d u0 u1 u2 p', **{'cls': GridVariable})
# Find the turbule u, v, w
initial_equations += [Eq(u, ubar + amp * (lx / 2.0) * cx * sy * sz)]
initial_equations += [Eq(v, - amp * sx * cy * sz)]
initial_equations += [Eq(w, - amp * (lz / 2.0) * sx * sy * cz)]
initial_equations += [Eq(p, S.One / (gama * Minf**2.0))]
initial_equations += [Eq(d, S.One)]

# Now evaluate the Conservative vector

initial_equations += [Eq(q_vector[0], d)]
initial_equations += [Eq(q_vector[1], d * u), Eq(q_vector[2], d * v), Eq(q_vector[3], d * w)]
initial_equations += [Eq(q_vector[4], p / (gama - S.One) + 0.5 * d * (u**2 + v**2 + w**2))]

# Instantiate a grid based initialisation classes
initial = GridBasedInitialisation()
initial.add_equations(grid_equations + initial_equations)

# STEP 2
# Set the equation classes for the block (list)
block.set_equations([constituent, simulation_eq, initial] + stat_equation_classes)

# STEP 3
# Create the dictionary of schemes
schemes = {}
# Central scheme for spatial discretisation and add to the schemes dictionary
dsets = 'u0 u1 u2 T'
cent = StoreSome(4, dsets)
schemes[cent.name] = cent
# RungeKutta scheme for temporal discretisation and add to the schemes dictionary
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
# Set the discretisation schemes to be used (a python dictionary)
block.set_discretisation_schemes(schemes)

# STEP 4 add io for the block
kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays + [x, y, z]
output_hdf5 = iohdf5(save_every=100000, arrays=output_arrays, **kwargs)
block.setio([output_hdf5])

# STEP 6
# Perform the symbolic discretisation of the equations
block.discretise()

# STEP 7
# create an algorithm from the numerical solution
# Add a diagnostics class to monitor the simulation
arrays = ['rho', 'u0', 'p', 'T', 'u1']
arrays = [block.location_dataset('%s' % dset) for dset in arrays]
indices = [(0, 10, 20), (30, 40, 50), (50, 129, 50), (64, 'block0np0-1', 64), (64, 96, 64)]
SM = SimulationMonitor(arrays, indices, block, print_frequency=250, fp_precision=12, output_file='output.log')
alg = TraditionalAlgorithmRK(block, simulation_monitor=SM)

# STEP 8
# set the simulation data type: if not set "Double" is default
SimulationDataType.set_datatype(Double)

# STEP 9
# Write the OPSC compatible code for the numerical solution
OPSC(alg)

# STEP 10
# Populate the values of the constants like Re, Pr etc and the number of points for the
# simulation etc. In the future reading thes from HDF5 would be provided
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1',
    'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', "c0", "c1", "c2", "lx0", "lx2", "Twall"]
values = ['180.0', '1.4', '0.01', '0.72', '0.00001', '100000', '256', '256', '256',
    '11.0/block0np0', '2.0/(block0np1-1)', '4.0/block0np2', '-1', '0', '0', "11.0", "4.0", "1.0"]
substitute_simulation_parameters(constants, values)

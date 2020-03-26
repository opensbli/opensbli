#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from sympy import sin, log, cos, pi, tanh
from opensbli.utilities.helperfunctions import substitute_simulation_parameters, print_iteration_ops

# STEP 0 Create the equations required for the numerical solution
# Problem dimension
ndim = 3
stats = True

sc1 = "**{\'scheme\':\'Teno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) - KD(_i,_j)*c_j )" % sc1
energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) - Dot(c_j, u_j))" % sc1
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
# Substitutions
substitutions = [stress_tensor, heat_flux]
constants = ["Re", "Pr", "gama", "Minf", "c_j"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**0.7))"
#viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

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

# Expand speed of sound add the expanded equations to the constituent relations
eqns = einstein_eq.expand(speed_of_sound, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Expand viscosity add the expanded equations to the constituent relations
eqns = einstein_eq.expand(viscosity, ndim, coordinate_symbol, substitutions, constants)
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

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

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

direction = 1
side = 0
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)
# Isothermal wall condition
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
#boundaries[direction][side] = IsothermalWallBC(direction, 0, wall_eqns)
boundaries += [IsothermalWallBC(direction, 0, wall_eqns)]
# Side 1 (top) boundary
side = 1
wall_eqns = [rhoE_wall]
boundaries += [IsothermalWallBC(direction, 1, wall_eqns)]


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
#grid_equations += [Eq(x, i * dx), Eq(y, -S.One + j * dy), Eq(z, k * dz)]
ny = ConstantObject('block0np1')
stretch = ConstantObject('stretch')
Ly = 2
stretched_eqn =  0.5*Ly*(1.0-((tanh(stretch*(1.0-2.0*(j/(ny-1.0)))))/(tanh(stretch))))-1.0
grid_equations += [Eq(x, i * dx), Eq(y, stretched_eqn), Eq(z, k * dz)]
#for eqn in grid_equations:
    #pprint(eqn)
#exit()
# Initial condition
initial_equations = []

lx, ly, lz = symbols('lx0:%s' % ndim, **{'cls': ConstantObject})
sx, sy, sz = symbols('sx0:%s' % ndim, **{'cls': GridVariable})
cx, cy, cz = symbols('cx0:%s' % ndim, **{'cls': GridVariable})
ubar, amp, vonkar, b = symbols('ubar amp vonkar b', **{'cls': GridVariable})

# Set the constants for initialisation
initial_equations += [Eq(b, 5.5), Eq(vonkar, 2.5)]

# Equation for umean
initial_equations += [Eq(ubar, Piecewise(((1 - Abs(y)) * Re, (1 - Abs(y)) * Re < 10.0),
                         (vonkar * log((1 - Abs(y)) * Re) + b, True)))]
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

# Transform equations for non uniform grid distribution
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (False, False)], 2)
simulation_eq.apply_metrics(metriceq)

# STEP 2
# Set the equation classes for the block (list)
block.set_equations([constituent, simulation_eq, initial, metriceq] + stat_equation_classes)

# STEP 3
# Create the dictionary of schemes
weno_order = 6
Avg = RoeAverage([0, 1])
LLF = LLFTeno(weno_order, averaging=Avg)
schemes = {}
schemes[LLF.name] = LLF
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKuttaLS(3, formulation='SSP')
#rk = RungeKuttaLS(3)
schemes[rk.name] = rk
block.set_discretisation_schemes(schemes)

# STEP 4 add io for the block
kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays + [x, y, z, DataObject('D11')]
output_hdf5 = iohdf5(save_every=50000, arrays=output_arrays, **kwargs)
block.setio([output_hdf5])

if stats: 

	statistics = [DataObject('u2u2mean'), DataObject('rhou2mean'), DataObject('rhou2u1mean'), DataObject('u0u0mean'), 
		DataObject('rhou1u0mean'), DataObject('E_mean'), DataObject('u1u0mean'), DataObject('u1u1mean'), DataObject('rhou2u0mean'),
		DataObject('rhou0mean'), DataObject('rhou1mean'), DataObject('pp_mean'), DataObject('rhou2u2mean'), DataObject('u2mean'), 
		DataObject('M_mean'), DataObject('u2u0mean'), DataObject('p_mean'), DataObject('a_mean'), DataObject('T_mean'), 
		DataObject('rhou0u0mean'), DataObject('rhomean'), DataObject('mu_mean'), DataObject('u2u1mean'), DataObject('TT_mean'), 
		DataObject('rhou1u1mean'), DataObject('u0mean'), DataObject('u1mean'), DataObject('D11')] 

        kwargs = {'iotype': "Write", 'name': "stats_TENO.h5"}
	stat_arrays = statistics
	stats_hdf5 = iohdf5(arrays=stat_arrays, **kwargs)
	block.setio([stats_hdf5])

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
constants = ['Re', 'gama', 'Minf', 'Pr', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1',
    'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', "c0", "c1", "c2", "lx0", "lx2", "stretch", 'harten', 'TENO_CT', 'eps']
values = ['190.71', '1.4', '0.0955', '0.7', '1.0', '0.0002', '250000', '129', '129', '129',
    '4.0*M_PI/block0np0', '2.0/(block0np1-1)', '(4.0*M_PI/3.0)/block0np2', '-1', '0', '0', "4.0*M_PI", "(4.0*M_PI/3.0)", "1.7", '0.25', '1e-7', '1e-15']
substitute_simulation_parameters(constants, values)
print_iteration_ops()

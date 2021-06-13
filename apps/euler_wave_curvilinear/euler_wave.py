#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import pi, sin, pprint

# Number of dimensions of the system to be solved
ndim = 2

sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

# define the compresible Navier-Stokes equations in Einstein notation
sc1 = "**{\'scheme\':\'Weno\'}"
# governing equations for NS
a = "Conservative(detJ * rho*U_j,xi_j,%s)" % sc1
mass = "Eq(Der(rho,t), - %s/detJ)" % (a)
a = "Conservative(detJ * (rhou_i*U_j + p*D_j_i), xi_j , %s)" % sc1
momentum = "Eq(Der(rhou_i,t) , -  %s/detJ)" % (a)
a = "Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)" % sc1
energy = "Eq(Der(rhoE,t), - %s/detJ)" % (a)

# substitutions
substitutions = []
constants_eqns = ["gama", "Minf"]
# formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE-rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"

einstein = EinsteinEquation()

# create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
einstein.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
# eqns = einstein.expand(metric_vel, ndim, coordinate_symbol, substitutions, constants)
# for eq in eqns:
#     einstein.optional_subs_dict[eq.lhs] = eq.rhs
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

# change coordinate symbol to curvilinear
coordinate_symbol = "xi"
base_eqns = [mass, momentum, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, metric_vel]
# Expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = einstein.expand(base, ndim, coordinate_symbol, substitutions, constants_eqns)
    if base==momentum:
        for no, b in enumerate(base_eqns[i]):
            base_eqns[i][no] = Eq(base_eqns[i][no].lhs, base_eqns[i][no].rhs)
    else:
        if base==energy:
            base_eqns[i] = Eq(base_eqns[i].lhs, base_eqns[i].rhs)

# expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = einstein.expand(CR, ndim, coordinate_symbol, substitutions, constants_eqns)

for eqn in base_eqns:
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)

# reset coordinate symbol
coordinate_symbol = "x"

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)
# Select the numerical schemes
weno_order = 5
Avg = RoeAverage([0, 1])
LLF = LLFWeno(weno_order, formulation='Z', averaging=Avg)
schemes = {}
schemes[LLF.name] = LLF
# cent = Central(4)
fns = 'u0 u1 T'
cent = StoreSome(4, fns)
schemes[cent.name] = cent
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
block.set_discretisation_schemes(schemes)

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
constants = ["gama"]
for con in constants:
    local_dict[con] = ConstantObject(con)
# Initial conditions
# Curvilinear wavy mesh
(i, j) = block.grid_indexes
(Sx, Sy) = 6.0, 6.0
(dx, dy) = block.deltas
Ax, Ay = 0.04, 0.04
omega = 0.25
Lx, Ly = 2.0, 2.0
x0 = Eq(DataObject('x0'), i*dx + Ax*sin(2.0*pi*omega)*sin(Sx*pi*j*dy/Ly), evaluate=False)
x1 = Eq(DataObject('x1'), j*dy + Ay*sin(2.0*pi*omega)*sin(Sy*pi*i*dx/Lx), evaluate=False)

pprint(x0)
pprint(x1)
d_in = parse_expr("Eq(GridVariable(d), 1+0.2*sin(pi*(DataObject(x0)+DataObject(x1))))", local_dict=local_dict)
u0_in = parse_expr("Eq(GridVariable(u0), 1.0)", local_dict=local_dict)
u1_in = parse_expr("Eq(GridVariable(u1), -0.5)", local_dict=local_dict)
p_in = parse_expr("Eq(GridVariable(p), 1.0)", local_dict=local_dict)


rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5* d *(u0**2+u1**2))", local_dict=local_dict)

initial = GridBasedInitialisation()
initial.add_equations([x0, x1, d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE])

boundaries = [[0, 0] for t in range(ndim)]
# Periodic boundaries on all sides
for direction in [0, 1]:
    for side in [0, 1]:
        boundaries[direction][side] = PeriodicBC(direction, side)

block.set_block_boundaries(boundaries)

# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=1000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1'),  DataObject('D00'), DataObject('D01'), DataObject('D10'), DataObject('D11')])
block.setio([h5])

# Set equations on the block and discretise
block.set_equations([constituent, initial, simulation_eq, metriceq])
block.discretise()
# Generate an algorithm and write out an OPS C code
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

# Add the simulation parameters to the C code
constants = ['gama', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0']
values = ['1.4', '0.0005', '5000', '128', '128', '2.0/(block0np0)', '2.0/(block0np1)']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho_B0', every=100)
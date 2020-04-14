#!/usr/bin/env python

# Import all the functions from opensbli
from opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters, print_iteration_ops
from sympy import pi, cos, sin

# Number of dimensions of the system to be solved
ndim = 2

#ndim = 3
sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

einstein_expasion = EinsteinEquation()

# Stress tensor and heat flux
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
# Substitutions
substitutions = [stress_tensor, heat_flux]
constants = ["Re", "Pr","gama", "Minf", "SuthT", "RefT"]

# Viscous momentum and energy components
visc_momentum = "Eq(Der(rhou_i, t), Der(tau_i_j,x_j))"
visc_momentum = einstein_expasion.expand(visc_momentum, ndim, coordinate_symbol, substitutions, constants)
visc_momentum = [metriceq.apply_transformation(v) for v in visc_momentum]

visc_energy = "Eq(Der(rhoE, t), -Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"
visc_energy = einstein_expasion.expand(visc_energy, ndim, coordinate_symbol, substitutions, constants)
visc_energy = metriceq.apply_transformation(visc_energy)

# Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs

a = "Conservative(detJ * rho*U_j,xi_j,%s)" % sc1
mass = "Eq(Der(rho,t), - %s/detJ)" % (a)
a = "Conservative(detJ * (rhou_i*U_j + p*D_j_i), xi_j , %s)" % sc1
momentum = "Eq(Der(rhou_i,t) , -  %s/detJ)" % (a)
a = "Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)" % sc1
energy = "Eq(Der(rhoE,t), - %s/detJ)" % (a)
# Substitutions
substitutions = []

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(0.7))"

einstein_expasion.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = einstein_expasion.expand(metric_vel, ndim, coordinate_symbol, substitutions, constants)
for eq in eqns:
    einstein_expasion.optional_subs_dict[eq.lhs] = eq.rhs

# change the symbol to xi as we will be using metrics
curvilinear_symbol = "xi"
simulation_eq = SimulationEquations()
eqns = einstein_expasion.expand(mass, ndim, curvilinear_symbol, substitutions, constants)

simulation_eq.add_equations(eqns)

eqns = einstein_expasion.expand(momentum, ndim, curvilinear_symbol, substitutions, constants)
for no, eq in enumerate(eqns):
	eqns[no] = Eq(eqns[no].lhs, eqns[no].rhs + visc_momentum[no].rhs)

simulation_eq.add_equations(eqns)

eqns = einstein_expasion.expand(energy, ndim, curvilinear_symbol, substitutions, constants)
eqns = Eq(eqns.lhs, eqns.rhs + visc_energy.rhs)
simulation_eq.add_equations(eqns)


constituent = ConstituentRelations()
eqns = einstein_expasion.expand(velocity, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = einstein_expasion.expand(pressure, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = einstein_expasion.expand(speed_of_sound, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = einstein_expasion.expand(temperature, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = einstein_expasion.expand(viscosity, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)

# # Metric velocity
# eqns = einstein_expasion.expand(metric_vel, ndim, curvilinear_symbol, substitutions, constants)
# constituent.add_equations(eqns)

# Create a simulation block
block = SimulationBlock(ndim, block_number=0)
# Select the numerical schemes
weno_order = 5
Avg = SimpleAverage([0, 1])
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


Rx, Ry, circ = 3.0, 6.0, (5.0*pi/12.0)

# Generate circular cylinder grid
xi, eta = block.deltas[0]*block.grid_indexes[0]/1.0, block.deltas[1]*block.grid_indexes[1]/1.0
x = OpenSBLIEq(DataObject('x0'), -1.0*(Rx-(Rx-1)*xi)*cos(circ*(2*eta-1.0)))
y = OpenSBLIEq(DataObject('x1'), 1.0*(Ry-(Ry-1)*xi)*sin(circ*(2*eta-1.0)))
grid_eqns = [x,y]

# parse the initial conditions
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial = GridBasedInitialisation()
initial.add_equations(grid_eqns + initial_equations)

boundaries = []
# Create boundaries, one for each side per dimension
q_vector = flatten(simulation_eq.time_advance_arrays)
boundaries = []
direction = 0
boundaries += [DirichletBC(direction, 0, initial_equations)]

# Symmetry sphere
boundaries += [SymmetryBC(direction, 1)]
# Far field boundary
direction, side = 1,0
boundaries += [ExtrapolationBC(direction, side, order=0)]
direction, side = 1,1
boundaries += [ExtrapolationBC(direction, side, order=0)]
# Set the boundaries for the block
block.set_block_boundaries(boundaries)


# set the IO class to write out arrays
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=100000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio([copy.deepcopy(h5)])

# Set equations on the block and discretise
block.set_equations([constituent, initial, simulation_eq, metriceq])
block.discretise()
# Generate an algorithm and write out an OPS C code
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

# Add the simulation parameters to the C code
constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'TENO_CT', 'eps']
values = ['100.0', '1.4', '3.0', '0.71', '0.001', '10000', '120', '160', '1.0/(block0np0-1)', '1.0/(block0np1-1)', '1.0e-6', '1.0e-8']
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho_B0')


#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from opensbli.core.weno_opensbli import *
from opensbli.core.teno import LLFTeno
from sympy import pi, tan
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters

ndim = 2

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

metriceq = MetricsEquation()
metriceq.genreate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

# Create an optional substitutions dictionary, this will be used to modify 
# the equations when parsed
optional_subs_dict = metriceq.metric_subs

sc1 = "**{\'scheme\':\'Teno\'}"
sc2 = "**{\'scheme\':\'Central\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t),-(theta*Conservative(detvel_j*rho,xi_j,%s)/detJ + (1-theta)*Skew(rho*detvel_j,xi_j,%s)/detJ))" % (sc1, sc2)
momentum = "Eq(Der(rhou_i,t),-(theta*Conservative(detJ * (rho*U_j*u_i + p*D_j_i), xi_j , %s)/detJ + (1-theta)*(Skew(detvel_j*rhou_i, xi_j, %s)/detJ + Der(detJ*p*D_j_i,xi_j, %s)/detJ)))" % (sc1, sc2, sc2)
energy = "Eq(Der(rhoE,t),-(theta*Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)/detJ + (1-theta)*(Skew(rhoE*detvel_j,xi_j,%s)/detJ + Conservative(detvel_j*p,xi_j,%s)/detJ)))" % (sc1, sc2, sc2)

# Substitutions
substitutions = []

# Define all the constants in the equations
constants = ["gama", "Minf"]

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"

einstein_expasion = EinsteinEquation()
einstein_expasion.optional_subs_dict = optional_subs_dict

detvel = "Eq(detvel_i, detJ*D_i_j*u_j)"
metric_vel = "Eq(U_i, D_i_j*u_j)"
# detvel2 = "Eq(detvel2_j, detvel_j*"
eqns = einstein_expasion.expand(metric_vel, ndim, coordinate_symbol, substitutions, constants)
eqns += einstein_expasion.expand(detvel, ndim, coordinate_symbol, substitutions, constants)
for eq in eqns:
    einstein_expasion.optional_subs_dict[eq.lhs] = eq.rhs

# change the symbol to xi as we will be using metrics
curvilinear_symbol = "xi"
simulation_eq = SimulationEquations()
eqns = einstein_expasion.expand(mass, ndim, curvilinear_symbol, substitutions, constants)

simulation_eq.add_equations(eqns)

eqns = einstein_expasion.expand(momentum, ndim, curvilinear_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = einstein_expasion.expand(energy, ndim, curvilinear_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

constituent = ConstituentRelations()
eqns = einstein_expasion.expand(velocity, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = einstein_expasion.expand(pressure, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = einstein_expasion.expand(speed_of_sound, ndim, curvilinear_symbol, substitutions, constants)
constituent.add_equations(eqns)

## Evaluation of Ducros sensor in 2D
cart = CoordinateObject('x_i')
cartesian_coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]

x0, x1 = cartesian_coordinates[0], cartesian_coordinates[1]
u0, u1 = DataObject('u0'), DataObject('u1')
epsilon = ConstantObject('epsilon')
theta = DataObject('theta')
vorticity = (CentralDerivative(u1, x0) - CentralDerivative(u0, x1))


divergence = "Eq(divergence, (KD(_i,_j)*Der(u_i, x_j)))"
divergence = einstein_expasion.expand(divergence, ndim, coordinate_symbol, substitutions, constants)

ducros = Eq(theta, divergence.rhs**2 / (divergence.rhs**2 + vorticity**2 + epsilon))
pprint(ducros)
ducros = metriceq.apply_transformation(ducros)
pprint(ducros)

constituent.add_equations(ducros)

# for eqn in constituent.equations:
#     pprint(eqn)
# exit()

latex = LatexWriter()
latex.open('./equations.tex', 'Simulation equations ')
simulation_eq.write_latex(latex)
latex.close()

# # for eqn in simulation_eq.equations:
# #     pprint(eqn)

# # print "\n\n\n"

# # for eqn in constituent.equations:
# #     pprint(eqn)
# exit()
# Do multiple orders chage here if required
schemes = {}
Avg = SimpleAverage([0, 1])
teno_order = 5
LLF = LLFTeno(teno_order, averaging=Avg)
schemes[LLF.name] = LLF
rk = RungeKutta(3)
schemes[rk.name] = rk
# For the discretisation of the metrics use central scheme
cent = Central(4)
schemes[cent.name] = cent

block = SimulationBlock(ndim, block_number=0)


local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

init_eq = []
coordinates = [DataObject('x%d' % i) for i in range(ndim)]
x, y = coordinates[0], coordinates[1]
dx, dy = block.deltas[0], block.deltas[1]
i,j = block.grid_indexes[0], block.grid_indexes[1]
npy = block.ranges[1][1]
# grid generation
x_local = GridVariable('x0')
lx, ly = ConstantObject('lx'), ConstantObject('ly')
local_dy = GridVariable('ldy')
# we take the length in x to be [0,l] with a flatplate for x/3.0 in [0 & 1) and (2,3]
# for x/3.0 in [1 and 2] there is a ramp of angle theta
shock_loc = 40.0
ramp_angle = -3.08
local_ly = [ExprCondPair(ly, x_local < shock_loc)]
local_ly += [ExprCondPair(ly + (x_local - shock_loc) * tan(ramp_angle*pi/180), True)]

# Create the equations , x local, x and y
init_eq += [Eq(x_local, dx * i)]
init_eq += [Eq(x, x_local)]
init_eq += [Eq(local_dy, (1/(npy - 1.0)) * Piecewise(*local_ly))]
init_eq += [Eq(y, local_dy*j)]

gama = ConstantObject('gama')
minf = ConstantObject('Minf')
# Simpler way to write the initialisation or equations
d,u,v,p = symbols('d, u0, u1, pr', **{'cls': GridVariable })

boundaries = [[0, 0] for t in range(ndim)]
d_in = Eq(d, 1.0)
u0_in = Eq(u, 1.0)
u1_in = Eq(v, 0.0)
p_in = Eq(p, 1.0/(gama*minf**2.0))

# The names in strings should match with the names from equations
r, ru, rv, re = symbols('rho, rhou0, rhou1, rhoE', **{'cls': DataObject})

rho = Eq(r, d)
rhou0 = Eq(ru, d*u)
rhou1 = Eq(rv, d*v)
rhoE = Eq(re, p/(gama-1) + 0.5* d *(u**2 + v**2))

init_eq += [d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE]
initial = GridBasedInitialisation()
initial.add_equations(init_eq)

# Left dirichlet
direction = 0
side = 0
inlet_eq = [d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBC(direction, side, inlet_eq)

# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = OutletTransferBC(direction, side)

# Bottom inviscid wall
direction = 1
side = 0
boundaries[direction][side] = SymmetryBC(direction, side)

# Top dirichlet shock condition
direction = 1
side = 1
boundaries[direction][side] = SymmetryBC(direction, side)

block.set_block_boundaries(boundaries)

kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=100000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio(copy.deepcopy(h5))

block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), copy.deepcopy(metriceq), initial])
block.set_discretisation_schemes(schemes)

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
constants = ['gama', 'Minf', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'lx', 'ly', 'TENO_CT', 'eps', 'epsilon']
values = ['1.4', '2.0', '0.1', '10', '400', '255', '300.0/(block0np0-1.0)', '115.0/(block0np1-1.0)', "300.0", "115.0", '1e-5', '1e-15', '1e-15']
substitute_simulation_parameters(constants, values)

    
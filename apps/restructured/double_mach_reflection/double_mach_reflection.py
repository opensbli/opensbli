#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from opensbli.core.weno_opensbli import *
import copy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from sympy import sin, cos, pi

ndim = 2

sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
a = "Conservative(rhou_j,x_j,%s)" % sc1
mass = "Eq(Der(rho,t), - %s)" % (a)
a = "Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s)" % sc1
momentum = "Eq(Der(rhou_i,t) , -%s  )" % (a)
a = "Conservative((p+rhoE)*u_j,x_j, %s)" % sc1
energy = "Eq(Der(rhoE,t), - %s  )" % (a)
# Substitutions
substitutions = []

# Define all the constants in the equations
constants = ["gama", "Minf"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"

simulation_eq = SimulationEquations()
eq = EinsteinEquation()
eqns = eq.expand(mass, ndim, coordinate_symbol, substitutions, constants)

simulation_eq.add_equations(eqns)

eqns = eq.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

eqns = eq.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

constituent = ConstituentRelations()
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

eqns = eq.expand(speed_of_sound, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

schemes = {}
weno_order = 5
# Averaging procedure to be used for the eigen system evaluation
Avg = RoeAverage([0, 1])
# LLF scheme
LLF = LLFWeno(weno_order, formulation='Z', averaging=Avg)
# Add to schemes
schemes[LLF.name] = LLF
rk = RungeKutta(3)
schemes[rk.name] = rk

block = SimulationBlock(ndim, block_number=0)

block.set_discretisation_schemes(schemes)

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
for con in constants:
    local_dict[con] = ConstantObject(con)
# Initial conditions
x0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
x1 = parse_expr("Eq(DataObject(x1), block.deltas[1]*block.grid_indexes[1])", local_dict=local_dict)

x_loc = Eq(GridVariable('x_loc'), 1.0/6.0 + DataObject('x1')/(3**0.5))
rho_in = Eq(GridVariable('rho_in'), Piecewise((8.0, DataObject('x0') < GridVariable('x_loc')), (1.4, True)))
u_in = Eq(GridVariable('u_in'), Piecewise((8.25*cos(pi/6.0), DataObject('x0') < GridVariable('x_loc')), (0, True)))
v_in = Eq(GridVariable('v_in'), Piecewise((-8.25*sin(pi/6.0), DataObject('x0') < GridVariable('x_loc')), (0, True)))
p_in = Eq(GridVariable('p_in'), Piecewise((116.5, DataObject('x0') < GridVariable('x_loc')), (1.0, True)))

rho = parse_expr("Eq(DataObject(rho), rho_in)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), rho_in*u_in)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), rho_in*v_in)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p_in/(gama-1) + 0.5* rho_in *(u_in**2+v_in**2))", local_dict=local_dict)

initial = GridBasedInitialisation()
initial.add_equations([x0, x1, x_loc, rho_in, u_in, v_in, p_in, rho, rhou0, rhou1, rhoE])

boundaries = [[0, 0] for t in range(ndim)]
# Left free-stream at x= 0, inlet conditions
# Inflow conditions
direction = 0
side = 0
boundaries[direction][side] = InletTransferBC(direction, side)
# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = OutletTransferBC(direction, side)
# Bottom inviscid wall
direction = 1
side = 0
boundaries[direction][side] = SymmetryBC(direction, side)
# Top time-dependent dirichlet shock condition
direction = 1
side = 1

dt = ConstantObject('dt')
current_iter = block.get_temporal_schemes[0].iteration_number
time_dependence = Eq(GridVariable('s'), 1.0/6.0 + (1+20*dt*current_iter)/(3**0.5))
rho_in = Eq(GridVariable('rho_in'), Piecewise((8.0, DataObject('x0') < GridVariable('s')), (1.4, True)))
u_in = Eq(GridVariable('u_in'), Piecewise((8.25*cos(pi/6.0), DataObject('x0') < GridVariable('s')), (0, True)))
v_in = Eq(GridVariable('v_in'), Piecewise((-8.25*sin(pi/6.0), DataObject('x0') < GridVariable('s')), (0, True)))
p_in = Eq(GridVariable('p_in'), Piecewise((116.5, DataObject('x0') < GridVariable('s')), (1.0, True)))

rho = parse_expr("Eq(DataObject(rho), rho_in)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), rho_in*u_in)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), rho_in*v_in)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p_in/(gama-1) + 0.5* rho_in *(u_in**2+v_in**2))", local_dict=local_dict)

moving_shock_eqns = [time_dependence, rho_in, u_in, v_in, p_in, rho, rhou0, rhou1, rhoE]
boundaries[direction][side] = DirichletBC(direction, side, moving_shock_eqns)

block.set_block_boundaries(boundaries)

kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')])
block.setio(copy.deepcopy(h5))

block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial])

block.discretise()

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
constants = ['gama', 'Minf', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0']
values = ['1.4', '10.0', '0.0002', '1000', '400', '400', '4.0/(block0np0-1)', '1.0/(block0np1-1)']
substitute_simulation_parameters(constants, values)

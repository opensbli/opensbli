#!/usr/bin/env python

# Inviscid Delery bump with forcing: Hiten Mulchandani 2019
# Import all the functions from opensbli
from opensbli import *
import copy
from sympy import atan, pi, sqrt, sin, sinh, And, exp, cos
import numpy
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from opensbli.utilities.katzer_init import Initialise_Katzer

#############################################################################################################################################
#																																			#
# Inputs			 																														#
#																																			#
#############################################################################################################################################

input_dict = {
    "gama"                 : "1.4",
    "Minf"                 : "0.85", # to get a shock over the bump
    "dt"                   : "0.01", # original number is 0.01
    "niter"                : "100000", # set as 1 to de-bug
    "block0np0"            : "701", # original number is 1401
    "block0np1"            : "1401", # original number is 1401
    "lx"                   : "700.0", # change this manually under domain boundaries (below) as well!
    "ly"                   : "700.0", # change this manually under domain boundaries (below) as well!
    "Delta0block0"         : "lx/(block0np0-1)",
    "Delta1block0"         : "ly/(block0np1-1)",
    "S"                    : "0.1772453851", # unsteady net mass source
    "sigma"                : "1.0", # parameter that controls the width of the Gaussian function
    "xs"                   : "500", # x-location of unsteady mass source
    "ys"                   : "50", # y-location of unsteady mass source
    "fs"                   : "0.007", # (non-dimensional) frequency of unsteady mass source
    "phi"                  : "0.0" # phase of unsteady mass source
}
constants = input_dict.keys()
values = input_dict.values()

#############################################################################################################################################
#																																			#
# Governing equations - homogeneous base equations 																							#
#																																			#
#############################################################################################################################################

# define number of dimensions
ndim = 2
# define coordinate direction symbol (x) this will be x_i, x_j
coordinate_symbol = "x"
# added metric class to genenerate curvilinear transformation
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

lx, ly = symbols("lx ly", **{'cls':ConstantObject})
CTD.add_constant([lx, ly])

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
pressure = "Eq(p, (gama-1)*(rhoE-rho*(1/2)*(KD(_i,_j)*u_i*u_j)))" # see SESA3029 L6.4 for proof
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)" # see workbook for proof
temperature = "Eq(T, gama*p*Minf*Minf/rho)" # see workbook for proof

#############################################################################################################################################
#                                                                                                                                           #
# Schemes                                                                                                                                   #
#                                                                                                                                           #
#############################################################################################################################################

schemes = {}
weno_order = 5
# averaging procedure to be used for the eigen system evaluation
Avg = SimpleAverage([0, 1])
# LLF scheme
LLF = LLFWeno(weno_order, formulation='Z', averaging=Avg)
# add to schemes
schemes[LLF.name] = LLF
central = Central(4) # central scheme for the metric terms
schemes[central.name] = central
rk = RungeKuttaLS(3)
schemes[rk.name] = rk

block = SimulationBlock(ndim, block_number=0)
block.set_discretisation_schemes(schemes)

#############################################################################################################################################
#                                                                                                                                           #
# Governing equations - forcing term(s)                                                                                                     #
#                                                                                                                                           #
#############################################################################################################################################

constituent = ConstituentRelations()
# acoustic forcing equations
S, sigma, xs, ys, fs, phi, dt = symbols('S sigma xs ys fs phi dt', **{'cls':ConstantObject})
x, y = [DataObject('x%d' % i) for i in range(ndim)]
current_iter = block.get_temporal_schemes[0].iteration_number
# source term is an unsteady (time-harmonic) mass displacement with a Gaussian term to position the source at (xs, ys)
forcing_eqn = Eq(DataObject('q'), (S/sqrt(pi*sigma))*exp(-(sqrt((x-xs)**2+(y-ys)**2)/sigma)**2)*cos(2*pi*fs*dt*current_iter+phi))
constituent.add_equations(forcing_eqn)
pprint(forcing_eqn)

#############################################################################################################################################
#                                                                                                                                           #
# Governing equations - add forcing term(s) to base equations                                                                               #
#                                                                                                                                           #
#############################################################################################################################################

einstein = EinsteinEquation()
# create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs
einstein.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
simulation_eq = SimulationEquations()
#constituent = ConstituentRelations()

# change coordinate symbol to curvilinear
coordinate_symbol = "xi"

base_eqns = [mass, momentum, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, metric_vel]
# expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = einstein.expand(base, ndim, coordinate_symbol, substitutions, constants_eqns)
    if base==mass:
        base_eqns[i] = Eq(base_eqns[i].lhs, base_eqns[i].rhs + DataObject('q'))
        pprint(base_eqns[i].rhs)
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

#############################################################################################################################################
#																																			#
# Grid and intial conditions																												#
#																																			#
#############################################################################################################################################

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
grid_const = ["lx", "ly"]

for i in range(len(grid_const)):
    local_dict[grid_const[i]] = ConstantObject(grid_const[i])
lx, ly = ConstantObject('lx'), ConstantObject('ly')
x0 = Eq(DataObject('x0'), block.deltas[0]*block.grid_indexes[0]) # for a flat plate and bump
x_local = x0.lhs # need this for Piecewise function that comes later on (which needs a LHS only)
j = block.grid_indexes[1]

# bump profile
# xbump values taken from A.A.Lawal's PhD thesis, pg. 165
xbump1 = Eq(GridVariable('xbump1'), (99.31/700.0)*lx)
xbump2 = Eq(GridVariable('xbump2'), (100.71/700.0)*lx)
xbump3 = Eq(GridVariable('xbump3'), (256.88/700.0)*lx)
xbump4 = Eq(GridVariable('xbump4'), (286.37/700.0)*lx)
xbump5 = Eq(GridVariable('xbump5'), (384.13/700.0)*lx)
xbump6 = Eq(GridVariable('xbump6'), (388.75/700.0)*lx)

xbump_eqns = [xbump1, xbump2, xbump3, xbump4, xbump5, xbump6]

xb1 = xbump1.lhs
xb2 = xbump2.lhs
xb3 = xbump3.lhs
xb4 = xbump4.lhs
xb5 = xbump5.lhs
xb6 = xbump6.lhs

# domain boundaries
lx = 700.0 # length of domain in x-direction
ly = 700.0 # length of domain in y-direction

np0, np1 = ConstantObject('block0np0'), ConstantObject('block0np1')

# define constants for Delery bump sections
conv = pi/180.0 # conversion from degrees to radians
theta1 = 270.0 # for part b
dtheta1 = 4.0/(xb2-xb1) # for part b
R1 = (20.0/700.0)*lx # for parts b and f
C1 = (10.92/700.0)*lx # for part c
C2 = (0.05/700.0)*lx # for part c
R2 = (500.0/700.0)*lx # for part d
C3 = (0.01/700.0)*lx # for parts d and e
C4 = (0.8692365/700.0)*lx # for part d
C5 = (1.0299/700.0)*lx # for part d
R3 = (422.67/700.0)*lx # for part e
max_height = (12.0/700.0)*lx # max height of bump, for parts d and e
theta2 = 256.6287 # for part f
dtheta2 = 13.3713/(xb6-xb5) # for part f

arg1 = Eq(GridVariable('arg1'), R1*sin((theta1-dtheta1*(x_local-xb1+C3))*conv)+R1)
ar1 = arg1.lhs
arg2 = Eq(GridVariable('arg2'), ((x_local-xb2)/(xb3-xb2))*C1+C2)
ar2 = arg2.lhs
arg3 = Eq(GridVariable('arg3'), ((sqrt(R2**2-((x_local-xb4+C3)**2))-R2)/C4)*C5+max_height)
ar3 = arg3.lhs
arg4 = Eq(GridVariable('arg4'), sqrt(R3**2-(x_local-xb4-C3)**2)-R3+max_height)
ar4 = arg4.lhs
arg5 = Eq(GridVariable('arg5'), R1*sin((theta2+dtheta2*(x_local-xb5+C3))*conv)+R1)
ar5 = arg5.lhs
arg_eqns = [arg1, arg2, arg3, arg4, arg5]

y0 = Eq(GridVariable('y0'), Piecewise((0.0, x_local<xb1), (ar1, And(x_local>=xb1, x_local<xb2)), (ar2, And(x_local>=xb2, x_local<xb3)), \
	(ar3, And(x_local>=xb3, x_local<xb4)), (ar4, And(x_local>=xb4, x_local<xb5)), (ar5, And(x_local>=xb5, x_local<xb6)), (0.0, True)))
ly0 = y0.lhs

local_y = Eq(DataObject('x1'), ly0+j*(ly-ly0)/(np1-1)) # for inviscid grid

# equations to evaluate for the grid
grid_equations = [x0] + xbump_eqns + arg_eqns + [y0, local_y]

# initial flow conditions
d_in = parse_expr("Eq(GridVariable(d), 1.0)", local_dict=local_dict)
u0_in = parse_expr("Eq(GridVariable(u0), 1.0)", local_dict=local_dict)
u1_in = parse_expr("Eq(GridVariable(u1), 0.0)", local_dict=local_dict)
p_in = parse_expr("Eq(GridVariable(p), 1.0/(gama*Minf**2.0))", local_dict=local_dict)  # p = 1/(gama*Minf^2)
rho = parse_expr("Eq(DataObject(rho), d)", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), d*u0)", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), d*u1)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), p/(gama-1) + 0.5*d*(u0**2+u1**2))", local_dict=local_dict)
flow_conditions = [d_in, u0_in, u1_in, p_in, rho, rhou0, rhou1, rhoE]

initial = GridBasedInitialisation()
initial.add_equations(grid_equations + flow_conditions)

#############################################################################################################################################
#																																			#
# Boundary conditions																														#
#																																			#
#############################################################################################################################################

boundaries = [[0, 0] for t in range(ndim)]

# left: inlet pressure extrapolation at inlet
direction = 0
side = 0
boundaries[direction][side] = InletPressureExtrapolateBC(direction, side)

# right: extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)

# bottom: inviscid wall
direction = 1
side = 0
boundaries[direction][side] = InviscidWall2DBC(direction, side) # for inviscid wall

# top: extrapolation boundary on top
direction = 1
side = 1
boundaries[direction][side] = ExtrapolationBC(direction, side, order=0)

# set boundary conditions
block.set_block_boundaries(boundaries)

#############################################################################################################################################
#																																			#
# Data i/o 																																	#
#																																			#
#############################################################################################################################################

# write data
save_every = 250
kwargs = {'iotype': "Write"}
#h5 = iohdf5(**kwargs)
h5 = iohdf5(save_every=save_every, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)
h5.add_arrays([DataObject('x0'), DataObject('x1')]) # save grid coordinates
h5.add_arrays([DataObject('p')]) # save pressure
#h5.add_arrays([DataObject('D00'), DataObject('D01'), DataObject('D10'), DataObject('D11')]) # save metric terms
#h5.add_arrays([DataObject('detJ')]) # save the determinant of the Jacobian matrix
block.setio(copy.deepcopy(h5))

#############################################################################################################################################
#																																			#
# Code generation																															#
#																																			#
#############################################################################################################################################

block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_eq), initial, metriceq])
block.discretise()
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)
substitute_simulation_parameters(constants, values)
print_iteration_ops(NaN_check='rho_B0')
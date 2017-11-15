#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from opensbli.core.weno_opensbli import *
from opensbli.core.teno import LLFTeno
from sympy import pi, tan
import copy
import re
import os
import shutil
from simulation_parameters2d import substitute_parameters

ndim = 2
# As the convective terms are to be modified for curvilinear transformation
# first generate the metrics

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

metriceq = MetricsEquation()
metriceq.genreate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

# Create an optional substitutions dictionary, this will be used to modify 
# the equations when parsed
optional_subs_dict = metriceq.metric_subs

shock_capturing_type = "weno"
sc1 = "**{\'scheme\':\'%s\'}" % shock_capturing_type.title()
# Define the compresible Navier-Stokes equations in Einstein notation.
a = "Conservative(detJ * rho*U_j,xi_j,%s)" % sc1
mass = "Eq(Der(rho,t), - %s/detJ)" % (a)
a = "Conservative(detJ * (rho*U_j*u_i + p*D_j_i), xi_j , %s)" % sc1
momentum = "Eq(Der(rhou_i,t) , -  %s/detJ)" % (a)
a = "Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)" % sc1
energy = "Eq(Der(rhoE,t), - %s/detJ)" % (a)
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
latex = LatexWriter()
latex.open('./equations.tex', 'Simulation equations ')
simulation_eq.write_latex(latex)
latex.close()
# Do multiple orders 
if shock_capturing_type.lower() == "weno":
    orders = [3, 5, '5Z']
elif shock_capturing_type.lower() == "teno":
    orders = [5, 6]

for order in orders:
    output_folder = "./%s_order_%s/" %(shock_capturing_type.lower(), str(order))
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)
    schemes = {}
    # Local LaxFredirich scheme for weno 
    # Averaging procedure to be used for the eigen system evaluation
    Avg = RoeAverage([0, 1])
    # LLF scheme
    # .title gives the first letter upercase
    #from sympy.core.compatibility import exec_
    exec('LLF = LLF%s(order, averaging=Avg)'%(shock_capturing_type.title()))
    # Add to schemes
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
    boundaries[direction][side] = DirichletBoundaryConditionBlock(direction, side, inlet_eq)

    # Right extrapolation at outlet
    direction = 0
    side = 1
    boundaries[direction][side] = OutletTransferBoundaryConditionBlock(direction, side)

    # Bottom inviscid wall
    direction = 1
    side = 0
    boundaries[direction][side] = SymmetryBoundaryConditionBlock(direction, side)

    # Top dirichlet shock condition
    direction = 1
    side = 1
    boundaries[direction][side] = SymmetryBoundaryConditionBlock(direction, side)

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
    constants = ['gama', 'Minf', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'lx', 'ly']
    values = ['1.4', '2.0', '0.1', '10000', '400', '255', '300.0/(block0np0-1.0)', '115.0/(block0np1-1.0)', "300.0", "115.0"]
    simulation_name = 'opensbli'
    if (shock_capturing_type.lower() == 'teno'):
        # TODO give the reference for these teno parameters
        # we need epsilon and TENO_CT values
        constants += ['eps', 'TENO_CT']
        values += ['1e-15']
        if order == 5:
            values += ['1e-5']
        elif order == 6:
            values += ['1e-7']
    substitute_parameters(simulation_name, constants, values)
    shutil.move('./latex_output', output_folder)
    for files in os.listdir('./'):
        if files.endswith((".cpp", ".h")):
            shutil.move(files, output_folder)

    
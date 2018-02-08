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
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
from opensbli.utilities.katzer_init import Initialise_Katzer

ndim = 2

# Type of shock capturing scheme weno or teno
shock_capturing_type = "weno"

# Type of Averaging *simple or roe
averaging_procedure = "simple"
# As the convective terms are to be modified for curvilinear transformation
# first generate the metrics

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

metriceq = MetricsEquation()
metriceq.genreate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

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


# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

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

latex = LatexWriter()
latex.open('./equations.tex', 'Simulation equations ')
simulation_eq.write_latex(latex)
latex.close()
# Do multiple orders chage here if required
if shock_capturing_type.lower() == "weno":
    orders = [5]
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
    exec('Avg = %sAverage([0, 1])' % (averaging_procedure.title()))
    # LLF scheme
    # .title gives the first letter upercase
    exec('LLF = LLF%s(order, averaging=Avg, formulation = "Z")'%(shock_capturing_type.title()))
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
    grid_const = ["lx", "ly"]
    for con in grid_const:
        local_dict[con] = ConstantObject(con)
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

    coordinate_eqns = init_eq[:]

    # Reynolds number, Mach number and free-stream temperature for the initial profile
    Re, xMach, Tinf = 950.0, 2.0, 288.0
    ## Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
    polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1'))]
    n_poly_coefficients = 50
    # grid_const = ["Lx1", "by"]
    for con in grid_const:
        local_dict[con] = ConstantObject(con)
    initial = Initialise_Katzer(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_eqns)

    # Left dirichlet
    boundaries = [[0, 0] for t in range(ndim)]
    direction = 0
    side = 0
    boundaries[direction][side] = InletTransferBC(direction, side)

    # Right extrapolation at outlet
    direction = 0
    side = 1
    boundaries[direction][side] = OutletTransferBC(direction, side)


    wall_const = ["Minf", "Twall"]
    for con in wall_const:
        local_dict[con] = ConstantObject(con)
    # Isothermal wall condition
    direction = 1
    side = 0
    rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
    wall_eqns = [rhoE_wall]
    boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns)

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
    constants = ['Twall', 'Re', 'Pr', 'gama', 'Minf', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'lx', 'ly', "SuthT", "RefT"]
    values = ['1.67619431', '950.0', '0.72', '1.4', '2.0', '0.1', '10000', '400', '255', '400.0/(block0np0-1.0)', '115.0/(block0np1-1.0)', "400.0", "115.0", '110.4', '288.0']
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
    substitute_simulation_parameters(constants, values)
    shutil.move('./latex_output', output_folder)
    for files in os.listdir('./'):
        if files.endswith((".cpp", ".h")):
            shutil.move(files, output_folder)

    
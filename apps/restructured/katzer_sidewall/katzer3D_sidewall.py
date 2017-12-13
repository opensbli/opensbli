#!/usr/bin/env python
from opensbli.core import *
from opensbli.utilities.katzer_init import *
from opensbli.core.weno_opensbli import SimpleAverage, RoeAverage, LLFWeno
from opensbli.core.teno import LLFTeno
from opensbli.physical_models.euler_eigensystem import *
from opensbli.initialisation import GridBasedInitialisation, iohdf5
from opensbli.utilities.helperfunctions import substitute_simulation_parameters
import copy

ndim = 3
sc1 = "**{\'scheme\':\'Weno\'}"
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s) + Der(tau_i_j,x_j) )" % sc1
energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s) - Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )" % sc1
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i) - (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (-mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
# Substitutions
substitutions = [stress_tensor, heat_flux]
constants = ["Re", "Pr", "gama", "Minf", "SuthT", "RefT"]
# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
speed_of_sound = "Eq(a, (gama*p/rho)**0.5)"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, (T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)))"

# Instatiate equation classes
eq = EinsteinEquation()
base_eqns = [mass, momentum, energy]
constituent_eqns = [velocity, pressure, speed_of_sound, temperature, viscosity]
# Expand the base equations
for i, base in enumerate(base_eqns):
    base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, substitutions, constants)
# Expand the constituent relations
for i, CR in enumerate(constituent_eqns):
    constituent_eqns[i] = eq.expand(CR, ndim, coordinate_symbol, substitutions, constants)

block = SimulationBlock(ndim, block_number=0)

local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

local_dict['Lx1'], local_dict['Lx2'] = ConstantObject('Lx1'), ConstantObject('Lx2')
local_dict['by'], local_dict['bz'] = ConstantObject('by'), ConstantObject('bz')

x_loc = parse_expr("Eq(GridVariable(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)

boundaries = [[0, 0] for t in range(ndim)]
# Left pressure extrapolation at x= 0, inlet conditions
direction = 0
side = 0
boundaries[direction][side] = InletTransferBC(direction, side)
# Right extrapolation at outlet
direction = 0
side = 1
boundaries[direction][side] = OutletTransferBC(direction, side)
# Bottom no-slip isothermal wall
direction = 1
side = 0
# boundaries[direction][side] = SymmetryBC(direction, side)
wall_const = ["Minf", "Twall"]
for con in wall_const:
    local_dict[con] = ConstantObject(con)
# Isothermal wall condition
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns)
# Top dirichlet shock generator condition
direction = 1
side = 1
rho = parse_expr("Eq(DataObject(rho), Piecewise((1.129734572, (x0)>40.0), (1.00000596004, True)))", local_dict=local_dict)
rhou0 = parse_expr("Eq(DataObject(rhou0), Piecewise((1.0921171, (x0)>40.0), (1.00000268202, True)))", local_dict=local_dict)
rhou1 = parse_expr("Eq(DataObject(rhou1), Piecewise((-0.058866065, (x0)>40.0), (0.00565001630205, True)))", local_dict=local_dict)
# rhou2 = parse_expr("Eq(DataObject(rhou2), 0.0)", local_dict=local_dict)
rhoE = parse_expr("Eq(DataObject(rhoE), Piecewise((1.0590824, (x0)>40.0), (0.94644428042, True)))", local_dict=local_dict)

dset = block.location_dataset('rhou2')
new_dset = increment_dataset(dset, 1, -1)
rhou2 = Eq(dset, new_dset)
pprint(rhou2.rhs.indices)
upper_eqns = [x_loc, rho, rhou0, rhou1, rhou2, rhoE]

# Split BC
# bcs = [DirichletBC(direction, side, upper_eqns), OutletTransferBC(direction, side)]
# boundaries[direction][side] = SplitBC(direction, side, bcs)

# Dirichlet only shock generator on top boundary
boundaries[direction][side] = DirichletBC(direction, side, upper_eqns)

direction = 2
side = 0
rhoE_wall = parse_expr("Eq(DataObject(rhoE), DataObject(rho)*Twall/(gama*(gama-1.0)*Minf**2.0))", local_dict=local_dict)
wall_eqns = [rhoE_wall]
boundaries[direction][side] = IsothermalWallBC(direction, side, wall_eqns)

direction = 2
side = 1
boundaries[direction][side] = SymmetryBC(direction, side)

# Create SimulationEquations and Constituent relations, add the expanded equations
simulation_eq = SimulationEquations()
constituent = ConstituentRelations()

for eqn in base_eqns:
    simulation_eq.add_equations(eqn)

for eqn in constituent_eqns:
    constituent.add_equations(eqn)
weno_order = 5
Avg = RoeAverage([0, 1])
LLF = LLFWeno(weno_order, formulation='Z', averaging=Avg)
schemes = {}
schemes[LLF.name] = LLF
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

metriceq = MetricsEquation()
# Stretched in wall normal and spanwise directions
metriceq.genreate_transformations(ndim, coordinate_symbol, [(False, False), (True, False), (True, False)], 2)
simulation_eq.apply_metrics(metriceq)

block.set_block_boundaries(boundaries)

sim_eq = copy.deepcopy(simulation_eq)
CR = copy.deepcopy(constituent)

# Perform initial condition
# Call the new polynomial based katzer initialisation, stretch factor 3 with 17 coefficients for the polynomial
Re, xMach, Tinf = 950.0, 2.0, 288.0
## Ensure the grid size passed to the initialisation routine matches the grid sizes used in the simulation parameters
polynomial_directions = [(False, DataObject('x0')), (True, DataObject('x1')), (True, DataObject('x2'))] # Stretched in x1, x2 directions, uniform in x0
n_poly_coefficients = 45
gridx0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
gridx1 = parse_expr("Eq(DataObject(x1), Lx1*sinh(by*block.deltas[1]*block.grid_indexes[1]/Lx1)/sinh(by))", local_dict=local_dict)
gridx2 = parse_expr("Eq(DataObject(x2), Lx2*sinh(bz*block.deltas[2]*block.grid_indexes[2]/Lx2)/sinh(bz))", local_dict=local_dict)
coordinate_evaluation = [gridx0, gridx1, gridx2]
initial = Initialise_Katzer(polynomial_directions, n_poly_coefficients, Re, xMach, Tinf, coordinate_evaluation)

# Arrays to write out to file
kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays
h5 = iohdf5(save_every=10000, arrays=output_arrays, **kwargs)

kwargs['name'] = 'grid.h5'
grid_h5 = iohdf5(arrays=[DataObject('x0'), DataObject('x1'), DataObject('x2')], **kwargs)
block.setio([h5, grid_h5])

# Set equations on the block and discretise
block.set_equations([sim_eq, CR, initial, metriceq])
block.set_discretisation_schemes(schemes)
block.discretise()


alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

constants = ['gama', 'Minf', 'Pr', 'Re', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2',
             'Delta0block0', 'Delta1block0', 'Delta2block0', 'SuthT', 'RefT', 'eps', 'TENO_CT', 'Lx1', 'Lx2', 'by', 'bz']
values = ['1.4', '2.0', '0.72', '950.0', '1.67619431', '0.05', '100', '400', '400', '100',
          '400.0/(block0np0-1)', '115.0/(block0np1-1)', '57.5/(block0np2-1)', '110.4', '288.0', '1e-15', '1e-5', '115.0', '57.5', '3.0', '3.0']
substitute_simulation_parameters(constants, values)

#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
import copy


# Problem dimension
ndim = 1

# Define the wave equation in Einstein notation.
wave = "Eq(Der(phi,t), -c_j*Der(phi,x_j))"

equations = [wave]

# Substitutions
substitutions = []

# The wave speed
constants = ["c_j"]

# Coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

simulation_eq = SimulationEquations()
eq = EinsteinEquation()
eqns = eq.expand(wave, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

constituent = ConstituentRelations()

block = SimulationBlock(ndim, block_number=0)
block.sbli_rhs_discretisation = True

boundaries = []
# Create boundaries, one for each side per dimension
for direction in [0]:
    boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
    boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}

block.set_block_boundaries(boundaries)

# Initial conditions
local_dict = {"block": block, "GridVariable": GridVariable, "DataObject": DataObject}
x0 = parse_expr("Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])", local_dict=local_dict)
phi = parse_expr("Eq(DataObject(phi), sin(2.0*pi*DataObject(x0)))", local_dict=local_dict)
initial = GridBasedInitialisation()
initial.add_equations([x0, phi])

kwargs = {'iotype': "Write"}
output_arrays = simulation_eq.time_advance_arrays + [DataObject('x0')]
h5 = iohdf5(arrays=output_arrays, **kwargs)

simulation = copy.deepcopy(simulation_eq)
block.set_equations([initial, simulation])
block.setio(copy.deepcopy(h5))

schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

block.set_discretisation_schemes(schemes)
# Discretise the equations on block
block.discretise()

# Algorithm for the block
alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)

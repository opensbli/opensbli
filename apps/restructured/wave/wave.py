#!/usr/bin/env python
import sys
from math import ceil

# Import local utility functions
from sympy import *
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *
from opensbli.initialisation import *

def dt(dx, c):
    """ Given a grid spacing dx and the wave speed c, return the value of dt such that the CFL condition is respected. """
    courant_number = 0.2
    return (dx*courant_number)/c

BUILD_DIR = os.getcwd()

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
eq = Equation()
eqns = eq.expand(wave, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

constituent = ConstituentRelations()

block= SimulationBlock(ndim, block_number = 0)
block.sbli_rhs_discretisation = True

boundaries = []
# Create boundaries, one for each side per dimension
for direction in range(ndim):
    boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
    boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]

block.set_block_boundaries(boundaries)

# Initial conditions
local_dict = {"block" : block, "GridVariable" : GridVariable, "DataObject" : DataObject}
phi = parse_expr("Eq(DataObject(phi), sin(2.0*pi*(block.grid_indexes[0])*block.deltas[0]))", local_dict = local_dict)
initial = GridBasedInitialisation()
initial.add_equations([phi])

kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every =100, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays)

simulation = copy.deepcopy(simulation_eq)
block.set_equations([simulation, initial])
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

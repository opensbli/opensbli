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

metriceq = MetricsEquation()
metriceq.genreate_transformations(ndim, coordinate_symbol, [(True, False)], 1)


simulation_eq = SimulationEquations()
eq = Equation()
eqns = eq.expand(wave, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
simulation_eq.apply_metrics(metriceq)

constituent = ConstituentRelations()
# Output equations in LaTeX format.
latex = LatexWriter()
latex.open(path=BUILD_DIR + "/equations.tex")
metadata = {"title": "Equations", "author": "", "institution": ""}
latex.write_header(metadata)
temp = flatten(simulation_eq.equations)
latex.write_expression(temp)
latex.write_footer()
latex.close()

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
phi = parse_expr("Eq(DataObject(phi), sin(2*pi*(block.grid_indexes[0])*block.deltas[0]))", local_dict = local_dict)
initial = GridBasedInitialisation()
initial.add_equations([phi])


simulation = copy.deepcopy(simulation_eq); CR = copy.deepcopy(constituent); metric = copy.deepcopy(metriceq)
block.set_equations([simulation , CR, metric, initial])


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

#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *


BUILD_DIR = os.getcwd()

# Problem dimension
ndim = 1

# Define the wave equation in Einstein notation.
wave = "Eq(Der(phi,t), -c_j*Der(phi,x_j, **{\'scheme\':\'Weno\'}))"
# wave = "Eq(Der(phi,t), -c_j*Der(phi,x_j, **{\'scheme\':\'Central\'}))"
# wave = "Eq( Der(phi,t) + Der(rho,t), -c_j*Conservative(phi*phi,x_j, **{\'scheme\':\'Weno\'}) - c_j*Conservative(phi*phi*phi, x_j, **{\'scheme\':\'Weno\'}))"

mass = []

# equations = [burg]

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
#pprint(simulation_eq.equations)
# Discretise the constituent relations if it contains any derivatives
#simulation_eq.add_constituent_relations(constituent)

block= SimulationBlock(ndim, block_number = 0)

## Create eigensystem
Euler = EulerEquations(ndim)
ev_dict, LEV_dict, REV_dict = Euler.generate_eig_system()
weno_order = 3
flat_eqns = flatten(simulation_eq.equations)
# GLF = GLFCharacteristic(ev_dict, LEV_dict, REV_dict, weno_order)

cart = CoordinateObject('%s_i'%(coordinate_symbol))
coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
speeds = [ConstantObject('c%d'%d) for d in range(ndim)]
speeds_dict = dict(zip(coordinates, speeds))



schemes = {}
SF = ScalarLocalLFScheme(weno_order, flat_eqns, speeds_dict, ndim)
schemes[SF.name] = SF
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

deriv = SF.grouped_eqns[0][0]
fn = deriv.args[0]
direction = deriv.args[-1]
direction_index = coordinates.index(direction)
pprint(direction_index)
pprint(fn)
pprint(direction)
side = -1
# SF.generate_weno(fn, direction, side)
interpolations = SF.pre_process(direction, direction_index)

SF.update_WenoSolutionType(interpolations)

pprint(interpolations)
local_eval, reconstruction = SF.post_process(interpolations)



exit()
block.sbli_rhs_discretisation = True
boundaries = [PeriodicBoundaryConditionBlock()]*2*ndim
block.set_block_boundaries(boundaries)
block.set_equations([simulation_eq])
block.set_discretisation_schemes(schemes)
block.discretise()

TraditionalAlgorithm(block)
# Create a descritisation class where the equations are deepcopied using deepcopy
#ex = copy.deepcopy(simulation_eq)
# block= SimulationBlock(ndim, block_number = 0)





#################

latex = LatexWriter()
latex.open('./equations.tex')
metadata = {"title": "Einstein Expansion of equations", "author": "Jammy", "institution": ""}
latex.write_header(metadata)
for index, eq in enumerate(flatten(simulation_eq.equations)):
    if isinstance(eq, Equality):
        latex.write_expression(eq)
latex.write_footer()
latex.close()



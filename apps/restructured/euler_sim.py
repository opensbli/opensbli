#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.physical_models.euler_eigensystem import *


BUILD_DIR = os.getcwd()

ndim = 2

weno = True
Euler_eq = EulerEquations(ndim, weno)

substitutions = []


coordinate_symbol = "x"

simulation_eq = SimulationEquations()
eq = Equation()

#WARNING: Not sure if the pressure in the momentum equation is being correctly evaluated 
# with the kronecker delta object, there is no pressure in the second parsed equation
for equation in Euler_eq.equations:
	eqns = eq.expand(equation, ndim, coordinate_symbol, substitutions, constants=Euler_eq.constants)
	pprint(eqns)
	print "-----------------------------------------------------------------"
	simulation_eq.add_equations(eqns)
constituent = ConstituentRelations()

for no, equation in enumerate(Euler_eq.formulas):
	pprint(no)
	print "Equation is: %s " % equation
	eqns = eq.expand(equation, ndim, coordinate_symbol, substitutions, constants=Euler_eq.constants)
	constituent.add_equations(eqns)
	print "------------------------------------------------------------------"
pprint(constituent.equations)

block= SimulationBlock(ndim, block_number = 0)

## Create eigensystem
ev_dict, LEV_dict, REV_dict = Euler_eq.generate_eig_system()
weno_order = 3
flat_eqns = flatten(simulation_eq.equations)
print "printing flat equations"
for item in flat_eqns:
	pprint(item)
exit()

cart = CoordinateObject('%s_i'%(coordinate_symbol))
coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
speeds = [ConstantObject('c%d'%d) for d in range(ndim)]
speeds_dict = dict(zip(coordinates, speeds))



schemes = {}
GLF = GLFCharacteristic(ev_dict, LEV_dict, REV_dict, weno_order)
schemes[GLF.name] = GLF
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

deriv = GLF.grouped_eqns[0][0]
fn = deriv.args[0]
direction = deriv.args[-1]
direction_index = coordinates.index(direction)
# pprint(direction_index)
# pprint(fn)
# pprint(direction)
side = -1
# SF.generate_weno(fn, direction, side)
interpolations = SF.pre_process(direction, direction_index)

SF.update_WenoSolutionType(interpolations)

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



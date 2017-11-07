#!/usr/bin/env python

""" This setup generates code required for the implementation of the error indicator described here:

C. T. Jacobs, M. Zauner, N. De Tullio, Satya P. Jammy, David J. Lusher, N. D. Sandham (Submitted). "An error indicator for finite difference methods using spectral techniques with application to aerofoil simulation". Submitted to the ParCFD 2017 Special Issue of Computers & Fluids.

and its application to a V2C airfoil simulation. """

import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.initialisation import *


BUILD_DIR = os.getcwd()

# Problem dimension
ndim = 2
coordinate_symbol = "x"

# Block definition.
block = SimulationBlock(ndim, block_number = 0)
local_dict = {"block" : block, "GridVariable" : GridVariable, "DataObject" : DataObject}

# Compute metrics.
metric_equation = MetricsEquation()
metric_equation.genreate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 1)
transformed = metric_equation.classical_strong_differentiabilty_transformation[0]
transformed = transformed.subs(metric_equation.general_function, DataObject("u0"))
pprint(transformed)

# Simulation equations.
eq = Equation()
vorticity = "Eq(Der(vorticity_i, t), LC(_i,_j,_k)*Der(velocity_k, x_j))"
simulation_equations = SimulationEquations()
eqns = eq.expand(vorticity, ndim, coordinate_symbol, substitutions=[], constants=[])
simulation_equations.add_equations(eqns)

# Apply metrics to simulation equations.
simulation_equations.apply_metrics(metric_equation)

constituent = ConstituentRelations()
velocity = "Eq(vort_i, vorticity_i)"
eqns = eq.expand(velocity, ndim, coordinate_symbol, substitutions=[], constants=[])
constituent.add_equations(eqns)

# Compute velocity from rhou_i/rho
initial = GridBasedInitialisation()
x0 = "Eq(DataObject(x0), block.deltas[0]*block.grid_indexes[0])"
x1 = "Eq(DataObject(x1), block.deltas[1]*block.grid_indexes[1])"
velocity0 = "Eq(DataObject(velocity0), DataObject(rhou0)/DataObject(rho))"
velocity1 = "Eq(DataObject(velocity1), DataObject(rhou1)/DataObject(rho))"
vorticity0 = "Eq(DataObject(vorticity0), 0)"
vorticity1 = "Eq(DataObject(vorticity1), 0)"
eqns = [velocity0, velocity1, vorticity0, vorticity1]
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
initial.add_equations(initial_equations)

# Numerical schemes.
schemes = {}
cent = Central(4)
schemes[cent.name] = cent
rk = RungeKutta(3)
schemes[rk.name] = rk

# Boundary conditions.
boundaries = []
direction = 0
boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]

direction = 1
boundaries += [PeriodicBoundaryConditionBlock(direction, 0)]
boundaries += [PeriodicBoundaryConditionBlock(direction, 1)]

block.set_block_boundaries(boundaries)
metric = copy.deepcopy(metric_equation)
block.set_equations([copy.deepcopy(constituent), copy.deepcopy(simulation_equations), metric, initial])
block.set_discretisation_schemes(schemes)
block.discretise()

# Quantities to initialise from HDF5 file.
rho = "Eq(DataObject(rho), r)"
rhou0 = "Eq(DataObject(rhou0), r*u0)"
rhou1 = "Eq(DataObject(rhou1), r*u1)"
rhoE = "Eq(DataObject(rhoE), p/(gama-1) + 0.5* r*(u0**2 + u1**2))"
eqns = [rho, rhou0, rhou1, rhoE]
initial_equations = [parse_expr(eq, local_dict=local_dict) for eq in eqns]
pprint(initial_equations)
arrays = [eq.lhs for eq in initial_equations + [parse_expr(x0, local_dict=local_dict), parse_expr(x1, local_dict=local_dict)]]
new_arays = block.dataobjects_to_datasets_on_block(arrays)
print block.block_datasets
for ar in new_arays :
    print str(ar)
    if "%s"%(ar) in block.block_datasets.keys():
        block.block_datasets["%s"%(ar)].read_from_hdf5 = True
        print         type(block.block_datasets["%s"%(ar)])
    else:
        block.block_datasets["%s"%(ar)] = ar
        block.block_datasets["%s"%(ar)].read_from_hdf5 = True
print block.block_datasets
#io = iohdf5(arrays=metric_equation.fdequations, iotype="read")

#pprint(metric_equation.spatial_discretisation(schemes, block))
pprint(metric_equation.fdequations)

alg = TraditionalAlgorithmRK(block)
SimulationDataType.set_datatype(Double)
OPSC(alg)


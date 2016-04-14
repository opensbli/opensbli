#!/usr/bin/env python
import sys
from math import ceil

# Import local utility functions
import opensbli
from opensbli.problem import *
from opensbli.latex import LatexWriter
from opensbli.spatial import *
from opensbli.bcs import *
from opensbli.ics import *
from opensbli.grid import *
from opensbli.timestepping import *
from opensbli.io import *
from opensbli.opsc import *

def dt(dx, c):
    """ Given a grid spacing dx and the wave speed c, return the value of dt such that the CFL condition is respected. """
    courant_number = 0.2
    return (dx*courant_number)/c

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 1D wave propagation simulation...")
start_total = time.time()

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

# Metrics
metrics = [False]

# Formulas for the variables used in the equations
formulas = []

# Create the problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
expanded_equations = problem.get_expanded(problem.equations)
expanded_formulas = problem.get_expanded(problem.formulas)

# Output equations in LaTeX format.
latex = LatexWriter()
latex.open(path=BUILD_DIR + "/equations.tex")
metadata = {"title": "Equations", "author": "", "institution": ""}
latex.write_header(metadata)
temp = flatten(expanded_equations)
latex.write_expression(temp)
temp = flatten(expanded_formulas)
latex.write_expression(temp)
latex.write_footer()
latex.close()

# Discretise the equations
start = time.time()

spatial_scheme = Central(8) # Eighth-order central differencing in space.
temporal_scheme = RungeKutta(3) # Third-order Runge-Kutta time-stepping scheme.

# Create a numerical grid of solution points
length = [1.0]*ndim
np = [1000]*ndim
deltas = [length[i]/(np[i]) for i in range(len(length))]

grid = Grid(ndim,{'delta':deltas, 'number_of_points':np})

# Perform the spatial discretisation
spatial_discretisation = SpatialDiscretisation(expanded_equations, expanded_formulas, grid, spatial_scheme)

# Perform the temporal discretisation
constant_dt = True
temporal_discretisation = TemporalDiscretisation(temporal_scheme, grid, constant_dt, spatial_discretisation)

boundary = BoundaryClass(grid)

# Update periodic boundary conditions in all direction
for dim in range(ndim):
    p = periodicboundary()
    boundary = p.apply_boundary(boundary, grid, temporal_discretisation.prognostic_variables, dim)

# Initial conditions
initial_conditions = ["Eq(grid.work_array(phi), sin(2*M_PI*(grid.Idx[0])*grid.deltas[0]))"]
initial_conditions = GridBasedInitialisation(grid, initial_conditions)

# I/O save conservative variables at the end of simulation
io = FileIO(temporal_discretisation.prognostic_variables)

# Grid parameters like number of points, length in each direction, and delta in each direction
c0 = 0.5
deltat = dt(deltas[0], c0)
niter = ceil(1.0/deltat)
print "Iterations: %d" % niter
l1 = ['niter', 'c0', 'deltat', 'precision', 'name']
l2 = [niter, c0, deltat, "double", "wave"]

# Constants in the system
simulation_parameters = dict(zip(l1,l2))

# Generate the code.
opsc = OPSC(grid, spatial_discretisation, temporal_discretisation, boundary, initial_conditions, io, simulation_parameters)
#opsc.set_diagnostics_level(1)
#opsc.generate()
#opsc.translate()


end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

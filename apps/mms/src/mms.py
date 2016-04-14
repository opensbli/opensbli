#!/usr/bin/env python

""" This simulation uses the Method of Manufactured Solutions to test solution convergence.

An manufactured/'analytical' solution for phi is substituted into the 2D advection-diffusion equation.
The non-zero residual is then subtracted from the RHS via a source term, such that
the manufactured solution then becomes the actual solution of the modified equation.

For more information about the procedure, see:

Roache (2002). Code Verification by the Method of Manufactured Solutions. Journal of Fluids Engineering, 124(1), 4-10. doi: 10.1115/1.1436090

"""

import sys
from math import ceil

# Import local utility functions
import opensbli
from opensbli.problem import *
from opensbli.latex import LatexWriter
from opensbli.spatial import *
from opensbli.bcs import *
from opensbli.grid import *
from opensbli.timestepping import *
from opensbli.io import *
from opensbli.opsc import *

def dt(dx, velocity):
    """ Given a grid spacing dx and the velocity, return the value of dt such that the CFL condition is respected. """
    courant_number = 0.025
    return (dx*courant_number)/velocity

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 2D MMS simulation...")
start_total = time.time()

# Problem dimension
ndim = 2

# Define the advection-diffusion equation in Einstein notation.
advection_diffusion = "Eq( Der(phi,t), -Der(phi*u_j,x_j) + k*Der(Der(phi,x_j),x_j) - s )"

equations = [advection_diffusion]

# Substitutions
substitutions = []

# Define all the constants in the equations
constants = ["k", "u_j", "s"]

# Coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False]

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

spatial_scheme = Central(DEGREE) # Central differencing in space.
temporal_scheme = RungeKutta(3) # Third-order Runge-Kutta time-stepping scheme.

# Create a numerical grid of solution points
length = [2*pi]*ndim
np = [NUMBER_OF_POINTS]*ndim
deltas = [length[i]/np[i] for i in range(len(length))]
grid = Grid(ndim,{'delta':deltas, 'number_of_points':np})

# Insert the source term 's'. The analytical solution for phi is sin(x[0]).
temp = EinsteinTerm('x_j')
x = temp.get_array(temp.get_indexed(ndim))
source_value = [-cos(x[0])*cos(x[1]) - 1.5*cos(x[1])*sin(x[0]) - 0.5*sin(x[0])*sin(x[1])]
source = EinsteinTerm("s")
x_grid = dict(zip(x, [grid.Idx[0]*grid.deltas[0], grid.Idx[1]*grid.deltas[1]]))
source_value = [v.subs(x_grid) for v in source_value]
expanded_equations[0][0] = expanded_equations[0][0].subs(source, source_value[0])
print expanded_equations[0][0]

# Perform the spatial discretisation
spatial_discretisation = SpatialDiscretisation(expanded_equations, expanded_formulas, grid, spatial_scheme)

# Perform the temporal discretisation
constant_dt = True
temporal_discretisation = TemporalDiscretisation(temporal_scheme, grid, constant_dt, spatial_discretisation)

# Create Boundary class to store the boundary conditions on the grid
boundary = BoundaryClass(grid)

# Update periodic boundary conditions in all direction
for dim in range(ndim):
    p = periodicboundary()
    boundary = p.apply_boundary(boundary, grid, temporal_discretisation.prognostic_variables, dim)

# Initial conditions. Note that we can use x0 and x1 as defined below and start off with the manufactured solution as the initial condition, but we'll start off with a zero initial condition instead to make it more rigorous.
x0 = "(grid.Idx[0]*grid.deltas[0])"
x1 = "(grid.Idx[1]*grid.deltas[1])"
initial_conditions = ["Eq(grid.work_array(phi), 0)"]
initial_conditions = GridBasedInitialisation(grid, initial_conditions)

# I/O save conservative variables at the end of simulation
io = FileIO(temporal_discretisation.prognostic_variables)

# Grid parameters like number of points, length in each direction, and delta in each direction
deltat = dt(max(deltas), velocity=2.0) # NOTE: We'll use an over-estimate for the velocity here in case of over-shoots.
T = 100.0 # NOTE: Make sure that the simulation runs long enough to ensure a steady-state solution is reached.
niter = ceil(T/deltat)
print "Going to do %d iterations." % niter

u0 = 1.0
u1 = -0.5
k = 0.75
simulation_parameters = {"niter":niter, "k":k, "u0":u0, "u1":u1, "deltat":deltat, "precision":"double", "name":SIMULATION_NAME}

# Generate the code.
opsc = OPSC(grid, spatial_discretisation, temporal_discretisation, boundary, initial_conditions, io, simulation_parameters)

end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

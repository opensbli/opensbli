#!/usr/bin/env python
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
    courant_number = 0.05
    return (dx*courant_number)/velocity

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 2D advection-diffusion simulation...")
start_total = time.time()

# Problem dimension
ndim = 2

# Define the advection-diffusion equation in Einstein notation.
advection_diffusion = "Eq( Der(phi,t), -Der(phi*u_j,x_j) + k*Der(Der(phi,x_j),x_j) )" 

equations = [advection_diffusion]

# Substitutions
substitutions = []

# Define all the constants in the equations
constants = ["k", "u_j"]

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

spatial_scheme = Central(4) # Fourth-order central differencing in space.
temporal_scheme = RungeKutta(3) # Third-order Runge-Kutta time-stepping scheme.

# Create a numerical grid of solution points
length = [10.0]*ndim
np = [10]*ndim
deltas = [length[i]/np[i] for i in range(len(length))]
print "Grid spacing is: ", deltas
grid = Grid(ndim,{'delta':deltas, 'number_of_points':np})

# Perform the spatial discretisation
spatial_discretisation = SpatialDiscretisation(expanded_equations, expanded_formulas, grid, spatial_scheme)

# Perform the temporal discretisation
const_dt = True
temporal_discretisation = TemporalDiscretisation(temporal_scheme, grid, const_dt, spatial_discretisation)

# Apply Boundary conditions
bcs = [("periodic", "periodic"), ("periodic", "periodic")]
boundary = BoundaryConditions(bcs, grid, temporal_discretisation.prognostic_variables)

# Initial conditions.
x0 = "(grid.Idx[0]*grid.deltas[0])"
x1 = "(grid.Idx[1]*grid.deltas[1])"
initial_conditions = ["Eq(grid.work_array(phi), exp(-( pow({x}-{x0}, 2)/(2*pow({spread_x}, 2)) + pow({y}-{y0}, 2)/(2*pow({spread_y}, 2)) )))".format(x=x0, x0=5, y=x1, y0=5, spread_x=0.075, spread_y=0.075)]
initial_conditions = GridBasedInitialisation(grid, initial_conditions)

# I/O save conservative variables at the end of simulation
io = FileIO(temporal_discretisation.prognostic_variables)

# Grid parameters like number of points, length in each direction, and delta in each direction
deltat = dt(max(deltas), velocity=1.5) # NOTE: We'll use an over-estimate for the velocity here in case of over-shoots.
T = 50.0
niter = ceil(T/deltat)
print "Time-step size is %f" % deltat
print "Going to do %d iterations." % niter

u0 = 0.1
u1 = 0.0
k = 0.1
simulation_parameters = {"niter":niter, "k":k, "u0":u0, "u1":u1, "deltat":deltat, "precision":"double", "name":"gaussian_bump"}

# Generate the code.
opsc = OPSC(grid, spatial_discretisation, temporal_discretisation, boundary, initial_conditions, io, simulation_parameters)

end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

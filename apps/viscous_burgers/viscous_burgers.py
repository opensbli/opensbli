""" A viscous Burgers equation solver.
Simulation parameters are based on those used here: http://pauli.uni-muenster.de/tp/fileadmin/lehre/NumMethoden/WS0910/ScriptPDE/Burgers.pdf """

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

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 1D viscous Burgers simulation...")
start_total = time.time()

# Problem dimension
ndim = 1

# Define the viscous Burgers equation in Einstein notation. Note that the advection term is written in the so-called "conservation form".
viscous_burgers = "Eq( Der(phi,t), -c_j*Conservative(phi*phi,x_j) + d*Der(Der(phi,x_j),x_j) )"
equations = [viscous_burgers]

# Substitutions
substitutions = []

# Constants
constants = ["c_j", "d"]

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

spatial_scheme = Central(4) # Fourth-order central differencing in space.
temporal_scheme = RungeKutta(3) # Third-order Runge-Kutta time-stepping scheme.

# Create a numerical grid of solution points
length = [10.0]
np = [200]
deltas = [length[0]/np[0]]

grid = Grid(ndim,{'delta':deltas, 'number_of_points':np})

# Perform the spatial discretisation
spatial_discretisation = SpatialDiscretisation(expanded_equations, expanded_formulas, grid, spatial_scheme)

# Perform the temporal discretisation
constant_dt = True
temporal_discretisation = TemporalDiscretisation(temporal_scheme, grid, constant_dt, spatial_discretisation)

# Boundary condition
boundary_condition = PeriodicBoundaryCondition(grid)
for dim in range(ndim):
    # Apply the boundary condition in all directions.
    boundary_condition.apply(arrays=temporal_discretisation.prognostic_variables, boundary_direction=dim)

# Initial conditions
x = "(grid.Idx[0]*grid.deltas[0])"
initial_conditions = ["Eq(grid.work_array(phi), exp(-((%s-3)*(%s-3))))" % (x, x)]
initial_conditions = GridBasedInitialisation(grid, initial_conditions)

# I/O save conservative variables at the end of simulation
io = FileIO(temporal_discretisation.prognostic_variables)

# Grid parameters like number of points, length in each direction, and delta in each direction
deltat = 0.05
T = 1.8
niter = ceil(T/deltat)
c0 = 0.5 # 1/2
d = 0.02 # Diffusion coefficient.
l1 = ['niter', 'c0', 'd', 'deltat', 'precision', 'name']
l2 = [niter, c0, d, deltat, "double", "viscous_burgers"]

# Constants in the system
simulation_parameters = dict(zip(l1,l2))

# Generate the code.
opsc = OPSC(grid, spatial_discretisation, temporal_discretisation, boundary_condition, initial_conditions, io, simulation_parameters)

end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

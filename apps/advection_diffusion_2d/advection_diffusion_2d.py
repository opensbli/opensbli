#!/usr/bin/env python
import sys

# Import local utility functions
import autofd
from autofd.problem import *
from autofd.algorithm import *
from autofd.latex import LatexWriter
from autofd.system import *

autofd.LOG.info("Generating code for the 2D advection-diffusion simulation...")

# Problem dimension
ndim = 2

# Define the equations in Einstein notation.
equations = ["Eq(Der(c,t),- conser(c,x_j))"]


# Substitutions
substitutions = []

# Define all the constants in the equations
constants = ["Lx0", "Lx1", "niter"]

# Define coordinate direction symbol (x) this will be x_i, x_j,x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False]

# Formulas for the variables used in the equations
formulas = []

# Create the problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
expanded_equations, expanded_formulas = problem.expand()

# Prepare the algorithm
temporal_scheme = "RK" # Runge-Kutta time-stepping scheme
temporal_order = 3
constant_timestep = True
spatial_scheme = "central_diff"
spatial_order = 4
evaluation_count = 0
work_array = "wk"
bcs = [("periodic", "periodic")]*ndim # Boundary conditions. The tuple represents (left, right) ends, and there is one tuple per dimension.
language = "OPSC"
multiblock = False
explicit_filter = [False]*ndim

algorithm = Algorithm(temporal_scheme, temporal_order, constant_timestep, spatial_scheme, spatial_order, work_array, evaluation_count, bcs, language, multiblock, explicit_filter)

# Output equations in LaTeX format.
latex = LatexWriter()
latex.open(path=BUILD_DIR + "/equations.tex")
metadata = {"title": "Equations", "author": "Christian T. Jacobs", "institution": "University of Southampton"}
latex.write_header(metadata)
temp = flatten([e.expandedeq for e in expanded_equations])
latex.write_equations(temp)
latex.write_footer()
latex.close()

# Prepare the computational system and generate the code
start = time.time()
system = System()

Lx0 = "2.0"
Lx1 = "2.0"
simulation_parameters = {"name":"advection_diffusion_2d", "Lx0":Lx0, "Lx1":Lx1, "nx0p[blk]":"10", "nx1p[blk]":"10", "dt":"0.0001", "niter":"10000", "a1":["2.0/3.0", "5.0/12.0", "3.0/5.0"], "a2":["1.0/4.0", "3.0/20.0", "3.0/5.0"], "dx0":"%s/nx0p[blk]" % Lx0, "dx1":"%s/nx1p[blk]" % Lx1}

system.prepare(expanded_equations, expanded_formulas, algorithm, simulation_parameters)
end = time.time()
LOG.debug('The time taken to prepare the system in %d Dimensions is %.2f seconds.' % (problem.ndim, end - start))

# Translate and compile the generated code
LOG.info("Compiling generated code...")
system.compile(algorithm, simulation_parameters)

#!/usr/bin/env python
import sys

# Import local utility functions
import autofd
from autofd.problem import *
from autofd.latex import LatexWriter
from autofd.system import *

autofd.LOG.info("Generating code for the 2D Taylor-Green Vortex simulation...")
start_total = time.time()
# Problem dimension
ndim = 2

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t),- Conservative(rhou_j,x_j))"
momentum = "Eq(Der(rhou_i,t) ,-Conservative(rhou_i*u_j + p* KD(_i,_j),x_j) + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t),- Conservative((p+rhoE)*u_j,x_j) +Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
lev = "Eq(vort_i, (LC(_i,_j,_k)*Der(u_k,x_j)))"
test = "Eq(Der(phi,t),- c_j* Der(phi,x_j))"

equations = [mass, momentum, energy]
#equations = [test]

# Substitutions
stress_tensor = "Eq(tau_i_j, (mu)*(Der(u_i,x_j)+ Conservative(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j,  (mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr","gama","mu", "Minf", "C23", "c_j"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False]

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - (1/(2))*(u_j*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(2/3))"
formulas = [velocity, pressure, temperature, viscosity]
#formulas = []

# Create the TGV problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
expanded_equations, expanded_formulas = problem.expand()

# Output equations in LaTeX format.
latex = LatexWriter()
latex.open(path=BUILD_DIR + "/equations.tex")
metadata = {"title": "Equations", "author": "", "institution": ""}
latex_substitutions = {'gama':'\gamma', 'rhou':'\\rho u', 'rhoE':'\\rho E'}
latex.write_header(metadata)
temp = flatten([e.expanded for e in expanded_equations])
latex.write_expression(temp, substitutions=latex_substitutions)
temp = flatten([e.expanded for e in expanded_formulas])
latex.write_expression(temp, substitutions=latex_substitutions)
latex.write_footer()
latex.close()

# Discretise the equations
start = time.time()

scheme = "central"
order = 4
spatial_scheme = Scheme(scheme, order)
temporal_scheme = RungeKutta(3) # Third-order Runge-Kutta time-stepping scheme.

# Create a numerical grid of solution points
grid = Grid(ndim) # FIXME: A HDF5 file or a user input

# Perform the spatial discretisation
spatial_discretisation = SpatialDiscretisation(expanded_equations,expanded_formulas, grid, spatial_scheme)

# Perform the temporal discretisation
const_dt = True
temporal_discretisation = TemporalDiscretisation(temporal_scheme, grid, const_dt, spatial_discretisation)

# Apply Boundary conditions
bcs = [("periodic", "periodic"), ("periodic", "periodic")]
boundary = BoundaryConditions(bcs, grid, temporal_discretisation.conservative)

# Initial conditions
initial_conditions = ["Eq(grid.work_array(phi), sin((grid.Idx[0] - 1)*grid.deltas[0]))"]
initial_conditions = GridBasedInitialisation(grid, initial_conditions)

# I/O save conservative variables at the end of simulation
io = FileIO(temporal_discretisation.conservative)

# Grid parameters like number of points, length in each direction, and delta in each direction
length = [1.0]*ndim
np = [8]*ndim
deltas = [length[i]/np[i] for i in range(len(length)) ] # how to define them
nsteps = 100

# Constants in the system
simulation_parameters = {"name":"testing"}

# Generate the code.
GenerateCode(grid, spatial_discretisation, temporal_discretisation, boundary, initial_conditions, io, simulation_parameters)

end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

latex = LatexWriter()
latex.open(path=BUILD_DIR + "/computations.tex")
metadata = {"title": "Equations", "author": "", "institution": ""}
latex.write_header(metadata)
#temp = flatten([e.expanded for e in expanded_equations])
#latex.write_equations(temp)
#temp = flatten([e.expanded for e in expanded_formulas])
#latex.write_equations(temp)
latex.write_footer()
latex.close()

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

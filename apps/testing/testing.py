#!/usr/bin/env python
import sys

# Import local utility functions
import autofd
from autofd.problem import *
from autofd.algorithm import *
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
#formulas = [velocity, pressure, temperature]
formulas = []

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

# Solve the equations on a grid
#grid = numerical_grid()# A HDF5 file or a user input
start = time.time()

sch = "central"
order = 4
spatial_scheme = Scheme(sch,order)
temporal_scheme = Scheme("Forward", 1)
grid = NumericalGrid(tuple(symbols('nx0:%d' % ndim, integer = True)))
const_dt = True
# get the spatial solution
spatial_solution = SpatialSolution(expanded_equations,expanded_formulas, grid, spatial_scheme)
# get the temporal solution
temporal_soln = TemporalSolution(temporal_scheme, grid,const_dt,spatial_solution)
# Apply Boundary conditions
bcs = [("periodic", "periodic"), ("periodic", "periodic")]
boundary = BoundaryConditions(bcs, grid, temporal_soln.conservative)
# Initial conditions
initialize = ["Eq(grid.work_array(phi), sin((grid.Idx[0] -1)*grid.deltas[0]))"]
Ics = GridBasedInitialization(grid, initialize)
# I/O save conservative variables at the end of simulation
IO = Fileio(temporal_soln.conservative)
end = time.time()
LOG.debug('The time taken to prepare the system in %d Dimensions is %.2f seconds.' % (problem.ndim, end - start))

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
LOG.debug('The time taken for the entire process for %d Dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

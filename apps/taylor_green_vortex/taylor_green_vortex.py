#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
import opensbli
from opensbli.problem import *
from opensbli.latex import *
from opensbli.spatial import *
from opensbli.bcs import *
from opensbli.ics import *
from opensbli.grid import *
from opensbli.timestepping import *
from opensbli.io import *
from opensbli.diagnostics import Reduction
from opensbli.opsc import OPSC

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 3D Taylor-Green Vortex simulation...")
start_total = time.time()

# Problem dimension
ndim = 3

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j))"
momentum = "Eq(Der(rhou_i,t) , -Skew(rhou_i*u_j,x_j) - Der(p,x_i) + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
ke = "Eq(ke, rho*(1/2)*u_j*u_j)"
enstrophy = "Eq(enstrophy, (1/2)*rho*(LC(_i,_j,_k)*Der(u_k,x_j))**2)"
rhomean = "Eq(rhomean, rho)"
equations = [mass, momentum, energy]
diagnostics = [ke, enstrophy, rhomean]

# Substitutions
stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr","gama", "Minf", "c_j"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False]

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(u_j*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
formulas = [velocity, pressure, temperature]

# Create the TGV problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
expanded_equations = problem.get_expanded(problem.equations)
expanded_formulas = problem.get_expanded(problem.formulas)
expanded_diagnostics = problem.get_expanded(problem.expand(diagnostics))

# Output equations in LaTeX format.
latex = LatexWriter()
latex.open(path=BUILD_DIR + "/equations.tex")
metadata = {"title": "Equations", "author": "", "institution": ""}
latex.write_header(metadata)
temp = flatten(expanded_equations)
latex.write_expression(temp[0])
# NOTE: Due to a bug in the LaTeX writer, you must replace (1.0/Re) and (1.0/((gama-1)*Minf*Minf*Pr*Re)) by their floating-point values in order to print the momentum and energy equations.
#latex.write_expression(temp[1])
#latex.write_expression(temp[2])
temp = flatten(expanded_formulas)
latex.write_expression(temp)
temp = flatten(expanded_diagnostics)
latex.write_expression(temp)
latex.write_footer()
latex.close()

# Discretise the equations
start = time.time()

spatial_scheme = Central(4) # Fourth-order central differencing in space.
temporal_scheme = RungeKutta(3) # Third-order Runge-Kutta time-stepping scheme.

# Create a numerical grid of solution points
length = [2.0*pi*1.0]*ndim
np = [64]*ndim
deltas = [length[i]/np[i] for i in range(len(length)) ]
print(deltas)
print(np)
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
x = "Eq(grid.grid_variable(x), grid.Idx[0]*grid.deltas[0])"
y = "Eq(grid.grid_variable(y), grid.Idx[1]*grid.deltas[1])"
z = "Eq(grid.grid_variable(z), grid.Idx[2]*grid.deltas[2])"
u = "Eq(grid.grid_variable(u),sin(x)*cos(y)*cos(z))"
v = "Eq(grid.grid_variable(v),-cos(x)*sin(y)*cos(z))"
w = "Eq(grid.grid_variable(w), 0.0)"
p = "Eq(grid.grid_variable(p), 1.0/(gama*Minf*Minf)+ (1.0/16.0) * (cos(2.0*x)+cos(2.0*y))*(2.0 + cos(2.0*z)))"
r = "Eq(grid.grid_variable(r), gama*Minf*Minf*p)"
initial_conditions = [x,y,z,u,v,w,p,r,
                      "Eq(grid.work_array(rho), r)",
                      "Eq(grid.work_array(rhou0), r*u)" ,
                      "Eq(grid.work_array(rhou1), r*v)",
                      "Eq(grid.work_array(rhou2), 0.0)",
                      "Eq(grid.work_array(rhoE), p/(gama-1) + 0.5* r *(u**2+ v**2 + w**2))"]

initial_conditions = GridBasedInitialisation(grid, initial_conditions)
end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

# Diagnostics
start = time.time()
reduction_type = ["sum", "sum", "sum", "sum", "sum"] # List of reduction types. This should be the same length as expanded_diagnostics.
red_eq = []
red_eq.append([Reduction(grid, expanded_diagnostics, expanded_formulas, temporal_discretisation.prognostic_variables, spatial_scheme, reduction_type, 100)])

end = time.time()
LOG.debug('The time taken to prepare the reductions in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

# Fix the deltat for 64^3 grid size and scale for other grid sizes accordingly
deltat = 3.385*10**-3/ (np[0]/64)
niter = ceil(20/deltat)
save_every = ceil(niter/8)

# I/O save conservative variables at every save_every steps and at the end of simulation
io = FileIO(temporal_discretisation.prognostic_variables, save_every)

start = time.time()
# Grid parameters like number of points, length in each direction, and delta in each direction
l1 = ['niter', 'Re', 'Pr', 'gama', 'Minf', 'mu', 'precision', 'name', 'deltat']

print "Going to do %d iterations." % niter
l2 = [niter, 1600, 0.71, 1.4, 0.1, 1.0, "double", "taylor_green_vortex", deltat]
# Constants in the system
simulation_parameters = dict(zip(l1,l2))

# Generate the code.
OPSC(grid, spatial_discretisation, temporal_discretisation, boundary_condition, initial_conditions, io, simulation_parameters, red_eq)

end = time.time()
LOG.debug('The time taken to generate the code in  %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

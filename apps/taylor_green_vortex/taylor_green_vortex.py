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

def dt(dx, c):
    """ Given a grid spacing dx and the wave speed c, return the value of dt such that the CFL condition is respected. """
    courant_number = 0.2
    return (dx*courant_number)/c

BUILD_DIR = os.getcwd()

opensbli.LOG.info("Generating code for the 3D Taylor-Green vortex simulation...")
start_total = time.time()

# Problem dimension
ndim = 3

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t),- Conservative(rhou_j,x_j))"
momentum = "Eq(Der(rhou_i,t) ,-Conservative(rhou_i*u_j + p* KD(_i,_j),x_j) + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t),- Conservative((p+rhoE)*u_j,x_j) +Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
equations = [mass, momentum, energy]

# Substitutions
stress_tensor = "Eq(tau_i_j, mu*(Der(u_i,x_j)+ Conservative(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j,  k*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["gama", "mu", "k", "cv"]

# Coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False, False]

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE/rho - (1/2)*(u_j*u_j)))"
temperature = "Eq(T, p/((gama-1)*rho*cv))"
viscosity = "Eq(mu, T**(2/3))"
formulas = [velocity, pressure, temperature, viscosity]

# Create the problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
expanded_equations, expanded_formulas = problem.get_expanded()

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
L = 1
length = [2*pi*L]*ndim
np = [25]*ndim
deltas = [length[i]/np[i] for i in range(len(length))]

grid = Grid(ndim,{'delta':deltas, 'number_of_points':np})

# Perform the spatial discretisation
spatial_discretisation = SpatialDiscretisation(expanded_equations, expanded_formulas, grid, spatial_scheme)

# Perform the temporal discretisation
const_dt = True
temporal_discretisation = TemporalDiscretisation(temporal_scheme, grid, const_dt, spatial_discretisation)

# Apply Boundary conditions
bcs = [("periodic", "periodic"), ("periodic", "periodic"), ("periodic", "periodic")]
boundary = BoundaryConditions(bcs, grid, temporal_discretisation.prognostic_variables)

# Initial conditions
u0 = 1.0
p0 = 100.0
T0 = 300.0
x = "(grid.Idx[0]*grid.deltas[0])"
y = "(grid.Idx[1]*grid.deltas[1])"
z = "(grid.Idx[2]*grid.deltas[2])"
rhou0 = "(%f*sin(%s)*cos(%s)*cos(%s))" % (u0, x, y, z)
rhou1 = "(-%f*cos(%s)*sin(%s)*cos(%s))" % (u0, x, y, z)
p = "%f + (1.0/16.0)*(cos(2*%s) + cos(2*%s))*(cos(2*%s) + 2.0)" % (p0, x, y, z)
ke = "0.5*(%s*%s + %s*%s)" % (rhou0, rhou0, rhou1, rhou1)
rhoE = "%s/(1.4-1.0) + %s" % (p, ke)

initial_conditions = ["Eq(grid.work_array(rho), 1.0)",
                      "Eq(grid.work_array(rhou0), %f*sin(%s)*cos(%s)*cos(%s))" % (u0, x, y, z),
                      "Eq(grid.work_array(rhou1), -%f*cos(%s)*sin(%s)*cos(%s))" % (u0, x, y, z),
                      "Eq(grid.work_array(rhou2), 0)",
                      "Eq(grid.work_array(rhoE), %s)" % rhoE,
                      "Eq(grid.work_array(p), %f + (1.0/16.0)*(cos(2*%s) + cos(2*%s))*(cos(2*%s) + 2.0))" % (p0, x, y, z),
                      "Eq(grid.work_array(T), %f)" % (T0)]
initial_conditions = GridBasedInitialisation(grid, initial_conditions)

# I/O save conservative variables at the end of simulation
io = FileIO(temporal_discretisation.prognostic_variables)

# Grid parameters like number of points, length in each direction, and delta in each direction
deltat = 3.385e-3
T = 5.0#10.0
niter = ceil(T/deltat)
#print "LOL", niter
mu = 1e-1
k = 0.05
cv = 718.0
gamma = 1.4
simulation_parameters = {"niter":niter, "gama":gamma, "mu":mu, "k":k, "cv":cv, "deltat":deltat, "precision":"double", "name":"taylor_green_vortex"}

# Generate the code.
opsc = OPSC(grid, spatial_discretisation, temporal_discretisation, boundary, initial_conditions, io, simulation_parameters)


end = time.time()
LOG.debug('The time taken to prepare the system in %d dimensions is %.2f seconds.' % (problem.ndim, end - start))

end_total = time.time()
LOG.debug('The time taken for the entire process for %d dimensions is %.2f seconds.' % (problem.ndim, end_total - start_total))

#!/usr/bin/env python
import sys

# Import local utility functions
import opensbli
from opensbli.problem import *
from opensbli.algorithm import *
from opensbli.latex import LatexWriter
from opensbli.system import *

opensbli.LOG.info("Generating code for the 2D Taylor-Green Vortex simulation...")

# Problem dimension
ndim = 2

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t),- Conservative(rhou_j,x_j))"
momentum = "Eq(Der(rhou_i,t) ,-Conservative(rhou_i*u_j + p* KroneckerDelta(_i,_j),x_j) + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t),-Conservative((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )"
equations = [mass, momentum, energy]

# Substitutions
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j) + Der(u_j,x_i) - (2/3) * KroneckerDelta(_i,_j)*Der(u_k,x_k)))"
heat_flux = "Eq(q_i, -(mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_i))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr", "mu", "gama", "Minf"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False]

# Formulas for the variables used in the equations

formulas = []

# Create the TGV problem and expand the equations.
problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
expanded_equations, expanded_formulas = problem.expand()

# Output equations in LaTeX format.
latex = LatexWriter()
latex.open(path=BUILD_DIR + "/equations.tex")
metadata = {"title": "Equations", "author": "Satya P Jammy", "institution": "University of Southampton"}
latex.write_header(metadata)
temp = flatten([e.expanded for e in expanded_equations])
latex.write_equations(temp)
latex.write_footer()
latex.close()

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

# Prepare the computational system and generate the code
start = time.time()
system = System()

Lx0 = "2.0*M_PI"
Lx1 = "2.0*M_PI"
simulation_parameters = {"name":"tgv_2d", "Lx0":Lx0, "Lx1":Lx1, "nx0p[blk]":"32", "nx1p[blk]":"32", "dt":"0.0005", "niter":"1000", "Minf":"0.1", "mu":"1.0", "a1":["2.0/3.0", "5.0/12.0", "3.0/5.0"], "a2":["1.0/4.0", "3.0/20.0", "3.0/5.0"], "gama":"1.4", "Pr":"0.71", "Re":"1600", "dx0":"%s/nx0p[blk]" % Lx0, "dx1":"%s/nx1p[blk]" % Lx1}

# Initial conditions
x = parse_expr('Eq(x, IDX0*dx0)', evaluate=False)
y = parse_expr('Eq(y, IDX1*dx1)', evaluate=False)
u_initial = parse_expr('Eq(u, sin(x)* cos(y))')
v_initial = parse_expr('Eq(v, -cos(x)* sin(y))')
rho_initial = parse_expr('Eq(rho, 1.0)', evaluate=False)
rhou0_initial = parse_expr('Eq(rhou0, sin(x)* cos(y))', evaluate=False)
rhou1_initial = parse_expr("Eq(rhou1, -cos(x)*sin(y))", evaluate=False)
rhoE_initial = parse_expr("Eq(rhoE, 1 + 0.5*(u**2 + v**2))", evaluate=False)            
            
initial_conditions = [x, y, u_initial, v_initial, rho_initial, rhou0_initial, rhou1_initial, rhoE_initial]

system.prepare(expanded_equations, expanded_formulas, algorithm, simulation_parameters, initial_conditions)
end = time.time()
LOG.debug('The time taken to prepare the system in %d Dimensions is %.2f seconds.' % (problem.ndim, end - start))

# Translate and compile the generated code
LOG.info("Compiling generated code...")
system.compile(algorithm, simulation_parameters)

#!/usr/bin/env python
import sys
import argparse

# Import local utility functions
import autofd
from autofd.problem import *
from autofd.algorithm import *
from autofd.latex import LatexWriter
from autofd.system import *

autofd.LOG.info("Generating code for the 2D Taylor-Green Vortex simulation...")

# Problem dimension
ndim = 2

# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t),- conser(rhou_j,x_j))"
momentum = "Eq(Der(rhou_i,t) ,-conser(rhou_i*u_j + p* KroneckerDelta(_i,_j),x_j) + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t),-conser((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )"
equations = [mass, momentum, energy]

# Substitutions
stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j) + Der(u_j,x_i) - (2/3) * KroneckerDelta(_i,_j)*Der(u_k,x_k)))"
heat_flux = "Eq(q_i, -(mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_i))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr", "mu", "gama", "Minf", "C23"]

# Define coordinate direction symbol (x) this will be x_i, x_j,x_k
coordinate_symbol = "x"

# Metrics
metrics = [False, False]

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - (1/(2*rho))*(rhou_j**2)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(2/3))"
formulas = [velocity, pressure, temperature, viscosity]

# Create the TGV problem and expand the equations.
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
#    latex = LatexWriter()
#    latex.open(path=BUILD_DIR + "/equations.tex")
#    metadata = {"title": "Equations", "author": "Satya P Jammy", "institution": "University of Southampton"}
#    latex.write_header(metadata)
#    temp = flatten([e.expandedeq for e in expanded_equations])
#    latex.write_equations(temp)
#    latex.write_footer()
#    latex.close()

# Prepare the computational system and generate the code
start = time.time()
system = System()
system.prepare(expanded_equations, expanded_formulas, algorithm)
end = time.time()
LOG.debug('The time taken to prepare the system in %d Dimensions is %.2f seconds.' % (problem.ndim, end - start))

# Translate and compile the generated code
LOG.info("Compiling generated code...")
system.compile(algorithm)

#!/usr/bin/env python
import sys

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

# Define coordinate direction symbol (x) this will be x_i, x_j,x_k
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


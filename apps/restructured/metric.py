#!/usr/bin/env python
import sys, os
from math import ceil

# Import local utility functions
#import opensbli as base
from opensbli.core import *
from opensbli.core.bcs import *
from opensbli.core.latex import *


BUILD_DIR = os.getcwd()

#base.LOG.info("Generating code for the 3D Taylor-Green Vortex simulation...")
#start_total = time.time()

# Problem dimension
ndim = 2
metric = {'scheme':'Central'}
# Define the compresible Navier-Stokes equations in Einstein notation.
mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,**%s))"%(metric)
momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j)  + Der(tau_i_j,x_j) )"
energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
ke = "Eq(ke, rho*(1/2)*Dot(u_j,u_j))"
enstrophy = "Eq(enstrophy, (1/2)*rho*(LC(_i,_j,_k)*Der(u_k,x_j))**2)"
rhomean = "Eq(rhomean, rho)"
equations = [ momentum, energy]
diagnostics = [ke, enstrophy, rhomean]

# Substitutions
stress_tensor = "Eq(tau_i_j, ((mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k))))"
heat_flux = "Eq(q_j, (1/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]

# Define all the constants in the equations
constants = ["Re", "Pr","gama", "Minf"]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(u_j*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
viscosity = "Eq(mu, T**(2/3))"
formulas = []
consta = {}

simulation_eq = MetricsEquation(ndim, coordinate_symbol, [(True, False), (False, False)], max_order = 2)
#pprint(simulation_eq.cart_to_curvilinear_functions)
#pprint(simulation_eq.curvilinear_to_cart_functions)
#sims = SimulationEquations()
eq = Equation()
coordinate_symbol = "xi"
#jjac = "Eq(jac_i_j, Der(x_i,xi_j))"
#second_ders = "Eq(sd_i_j_k, Der(jac_i_j,x_k))"
#metrics_eq = eq.expand(jjac , ndim, coordinate_symbol, substitutions, constants)
#pprint(metrics_eq)
#rhs = [m.rhs for m in metrics_eq]
#lhs = [m.lhs for m in metrics_eq]
#fdm = (Matrix(ndim,ndim, rhs))
#fds = Matrix(ndim,ndim, lhs)
#pprint(fdm)
#for no,fn in enumerate(simulation_eq.curvilinear_to_cart_functions):
    #if fn == 0:
        #fdm[no] = S.Zero
        #fds[no] = S.Zero
    #elif fn in simulation_eq.cartesian_coordinates:
        #fdm[no] = S.One
        #fds[no] = S.Zero
#pprint(fdm.adjugate())
#pprint(fdm.det())
#pprint(fds)
#subst = dict(zip(lhs, fds))
#pprint(subst)
#second_ders = "Eq(sd_i_j_k, Der(jac_i_j,xi_k))"
#second_eq = eq.expand(second_ders, ndim, coordinate_symbol, substitutions, constants)
#pprint(second_eq)
#second_eq = [s.subs(subst) for s in second_eq]
#pprint(second_eq)
metrics = [(True, True), (True, True)]
fd_fn = "Eq(fn_i_j, MetricDer(MetricDer(f,x_i), x_j))"
coordinate_symbol = "x"
fdfn = eq.expand(fd_fn, 2, coordinate_symbol, substitutions, constants, Metric=metrics)
original = "Eq(phi_i_m, "
#fds = Matrix([e.rhs for e in fdfn])
#pprint(fds)
latex = LatexWriter()
latex.open('./kernels3d.tex')
metadata = {"title": "Transformations", "author": "Jammy", "institution": ""}
latex.write_header(metadata)
latex.write_expression(fdfn)
latex.write_footer()
latex.close()
#latexify_expression(fdfn)
#pprint(fdfn)


#t = Transformation(ndim, coordinate_symbol, [(True, False), (True, False)])

#!/usr/bin/env python
# Import all the functions from opensbli
from opensbli import *
from sympy import sin, exp, pi, tan
import copy
from opensbli.postprocess  import *

ndim = 2
SimulationDataType.set_datatype(Double)

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"

metriceq = MetricsEquation()
metriceq.genreate_transformations(ndim, coordinate_symbol, [(True, True), (True, True)], 2)

# Velocity in 2D
vel = symbols("u0:%d"%ndim,  **{'cls':DataObject})
# Vorticity-z
wz = symbols("wz",  **{'cls':DataObject})
# coordinates
coord = symbols("x0:%d"%ndim,  **{'cls':CoordinateObject})

# Matrix of derivatives
der_matrix = Matrix(ndim,ndim,[CentralDerivative(u,x) for u in vel for x in coord])

vort = der_matrix[1,0] - der_matrix[0,1]
vort = metriceq.apply_transformation(vort)
vort = Eq(wz, vort)
post = PostProcessEquations()
post.add_equations(vort)

# Dilatation
dilatation = symbols("dl", **{'cls':DataObject})
dil = Eq(dilatation, der_matrix[0,0] + der_matrix[1,1])
post.add_equations(dil)

# We need only velocity
velocity = "Eq(u_i, rhou_i/rho)"
einstein_expasion = EinsteinEquation()
constituent = ConstituentRelations()
eqns = einstein_expasion.expand(velocity, ndim, coordinate_symbol, [], [])
constituent.add_equations(eqns)

latex = LatexWriter()
latex.open('./equations.tex', 'Postprocess equations ')
post.write_latex(latex)
latex.close()


block = SimulationBlock(ndim, block_number=0)


block.set_equations([post, constituent, metriceq])

schemes = {}
rk = RungeKutta(3)
schemes[rk.name] = rk
# For the discretisation of the metrics use central scheme
cent = Central(4)
schemes[cent.name] = cent


block.set_discretisation_schemes(schemes)


# BC's for block 
block0_bc = [[0, 0] for t in range(ndim)]

direction = 0
side = 0
block0_bc[direction][side] = DummyCarpenter(direction, side)

direction = 0
side = 1
block0_bc[direction][side] = DummyCarpenter(direction, side)


direction = 1
side = 0
block0_bc[direction][side] = DummyCarpenter(direction, side)

direction = 1
side = 1
block0_bc[direction][side] = DummyCarpenter(direction, side)

block.set_block_boundaries(block0_bc)


x,y = symbols("x0, x1", **{'cls':DataObject})
# set io
kwargs = {'iotype': "Write"}
h5 = iohdf5(**kwargs)
h5.add_arrays([wz, x, y, dilatation])
block.setio([h5])

kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([x, y] + list(symbols("rho rhou0 rhou1 rhoE",  **{'cls':DataObject})))
#wz = 
block.setio([h5_read])

block.discretise()
alg = TraditionalAlgorithmRK(block)
OPSC(alg)

# Write the simulation parameters

#!/usr/bin/env python
from opensbli import *
from sympy import sin, exp, pi, tan, cos
import copy
from opensbli.multiblock.algorithm import TraditionalAlgorithmRKMB
from sympy.functions.elementary.piecewise import Piecewise, ExprCondPair
import os
# Disable the gmpy library for this case to avoid deepcopy issues
os.environ['MPMATH_NOGMPY'] = '1'

import itertools
def create_exchange_calls_codes(multiblock_descriptor, dsets):
    kernels = []
    for block in multiblock_descriptor.blocks:
        arrays = [block.location_dataset(a) for a in dsets]
        kernels += block.apply_interface_bc(arrays, multiblock_descriptor)
    return kernels

ndim = 3
nblocks = 3
multi_block = MultiBlock(ndim, nblocks)
SimulationDataType.set_datatype(Double)

def generate_sponge_kernel(conserve_vector, block):
    """ Applies a sponge boundary on the outflow of the domain to damp oscillations."""
    length_sponge, lc, sigma = symbols("spongel, lc, sigma", **{'cls':GridVariable})
    gama, minf = symbols("gama Minf", **{'cls':ConstantObject})
    residual = symbols("Residual0:5", **{'cls':DataObject})
    x0 = symbols("x0", **{'cls':DataObject})
    equations = []
    equations += [Eq(length_sponge, 0.85)]
    # Characteristic length
    equations += [Eq(lc, x0 - (5.0 -length_sponge)), Eq(sigma,  0.5*(1.0 + cos(pi* lc/length_sponge)))]
    values = [1.0, 1.0, 0.0, 0.0, 1.0/(gama * minf**2.0 *(gama - 1.0)) + 0.5]
    for b0, b1, b2 in zip(residual, conserve_vector, values):
        equations += [Eq(b0, b0 - sigma * (b1- b2))]
    eqns = block.dataobjects_to_datasets_on_block(equations)
    # pprint ([eq for eq in eqns])
    ker = Kernel(block, computation_name="Sponge kernel block%d" %block.blocknumber)
    ker.kernelname = "sponge_kernel_block%d" %block.blocknumber
    ker.add_equation(eqns)
    ranges = copy.deepcopy(block.ranges)
    # for eqn in ker.equations:
    #     pprint(eqn)
    ker.ranges = ranges
    ker.ranges[0][0] =  ranges[0][1] - 62
    ker.update_block_datasets(block)
    return ker

def generate_wake_kernel(conserve_vector, mulitblock, wall_energy):
    """ Wake treatment at the block interface."""
    block = mulitblock.get_block(0)
    wk = symbols("wk0:5", **{'cls':DataObject})
    # Add the grid index if IDX ==0 then 
    # Also change the range of evaluation
    idx = block.grid_indexes[0]
    equations = [Eq(conserve_vector[0], 0.5 * (conserve_vector[0] + wk[0]))]
    # for rhou,v,w
    for b0, b1 in zip(conserve_vector[1:-1], wk[1:-1]):
        pairs = [ExprCondPair(0.0, Eq(idx,0)), ExprCondPair(0.5* (b0 + b1), True)]
        equations += [Eq(b0, Piecewise(*pairs, evaluate=False))]
    pairs = [ExprCondPair(wall_energy.rhs, Eq(idx,0)), ExprCondPair(0.5* (conserve_vector[-1] + wk[-1]), True)]
    equations += [Eq(conserve_vector[-1], Piecewise(*pairs, evaluate=False))]
    
    equations = block.dataobjects_to_datasets_on_block(equations)
    direction = 1
    side = 0
    # create it as a boundary condition, kernel, example we use DirichletBC
    bc = DirichletBC(direction, side, equations)
    ker = bc.apply([], block)
    ker.kernelname = "wake_treatment_kernel"
    ker.computation_name = "Wake treatment"
    ker.halo_ranges[1][0] = set()   
    # Wake exchanges from block2 wakeline (conserve_vector) to blokck0 work_arrays
    block2 = mulitblock.get_block(2)
    bc = InterfaceBC(direction, side,  match=(0, 1, 0, False))
    arrays = [block2.work_array(str(a)) for a in flatten(conserve_vector)]
    other_arrays = [block.work_array(str(a)) for a in flatten(wk)]
    wake_transfer1 = bc.apply_interface(arrays, block2, mulitblock, other_arrays=other_arrays)
    wake_transfer1.transfer_size[1] = 1
    wake_transfer1.transfer_from[1] = 0
    wake_transfer1.transfer_to[1] = 0
    wake_transfer1.computation_name = "waketransfer1"
    
    bc = InterfaceBC(direction, side,  match=(2, 1, 0, False))
    arrays = [block.work_array(str(a)) for a in flatten(conserve_vector)]
    wake_transfer2 = bc.apply_interface(arrays, block, mulitblock)
    wake_transfer2.transfer_size[1] = 1
    wake_transfer2.transfer_from[1] = 0
    wake_transfer2.transfer_to[1] = 0
    wake_transfer2.computation_name = "waketransfer2"
    return [wake_transfer1, ker, wake_transfer2]

# Define coordinate direction symbol (x) this will be x_i, x_j, x_k
coordinate_symbol = "x"
metriceq = MetricsEquation()
metriceq.generate_transformations(ndim, coordinate_symbol, [(True, True), (True, True), (False, False)], 2)
#Create an optional substitutions dictionary, this will be used to modify the equations when parsed
optional_subs_dict = metriceq.metric_subs

#Define the compresible Navier-Stokes equations in Einstein notation, by default the scheme is Central no need to
#Specify the schemes
mass = "Eq(Der(rho,t), - Skew(rho*u_j,x_j))"
momentum = "Eq(Der(rhou_i,t) , - Skew(rhou_i*u_j, x_j) - Der(p,x_i)  + Der(tau_i_j,x_j))"
energy = "Eq(Der(rhoE,t), - Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j) + Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j))"
# Substitutions used in the equations
stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
heat_flux = "Eq(q_j, (1.0/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
substitutions = [stress_tensor, heat_flux]
# Constants that are used
constants = ["Re", "Pr", "gama", "Minf", "mu"]
# Formulas for the variables used in the equations
velocity = "Eq(u_i, rhou_i/rho)"
pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
temperature = "Eq(T, p*gama*Minf*Minf/(rho))"

einstein_expasion = EinsteinEquation()
einstein_expasion.optional_subs_dict = optional_subs_dict

metric_vel = "Eq(U_i, D_i_j*u_j)"
eqns = einstein_expasion.expand(metric_vel, ndim, coordinate_symbol, substitutions, constants)
for eq in eqns:
    einstein_expasion.optional_subs_dict[eq.lhs] = eq.rhs

# Change the symbol to xi as we will be using metrics
coordinate_symbol = "x"
simulation_eq = SimulationEquations()

# Perform the expansion
eqns = einstein_expasion.expand(mass, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
eqns = einstein_expasion.expand(momentum, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)
eqns = einstein_expasion.expand(energy, ndim, coordinate_symbol, substitutions, constants)
simulation_eq.add_equations(eqns)

constituent = ConstituentRelations()
eqns = einstein_expasion.expand(velocity, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = einstein_expasion.expand(pressure, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)
eqns = einstein_expasion.expand(temperature, ndim, coordinate_symbol, substitutions, constants)
constituent.add_equations(eqns)

# Transform the equations into curvilinear form
simulation_eq.apply_metrics(metriceq)

# Set the equations on the blocks
multi_block.set_equations([simulation_eq, constituent, metriceq])

# Add filters to each block
filter_list = []
for no, block in enumerate(multi_block.blocks):
    if no == 1: # Main aerofoil block, C-mesh. Don't filter near the aerofoil
        j = block.grid_indexes[1]
        grid_condition = j >= 10
    else: # Filter everywhere 
        grid_condition = None
    filter_list += [BinomialFilter(block, order=10, grid_condition=grid_condition, sigma=0.01).equation_classes]
multi_block.set_filters(filter_list)

# Specify the numerical schemes
schemes = {}
rk = RungeKuttaLS(3)
schemes[rk.name] = rk
# cent = Central(4)
cent = StoreSome(4, 'u0 u1 u2 T')
schemes[cent.name] = cent
multi_block.set_discretisation_schemes(schemes)

# Initial conditions
d, u0, u1,u2, p = symbols("d, u0:3, p", **{'cls':GridVariable})
gama, Minf = symbols("gama, Minf", **{'cls':ConstantObject})
initial_equations = []
initial_equations += [Eq(d, 1.0)]
initial_equations += [Eq(u0, 1.0)]
initial_equations += [Eq(u1, 0.0)]
initial_equations += [Eq(u2, 0.0)]
initial_equations += [Eq(p, 1.0/(gama*Minf**2.0))]

# Set the conservative values
conserve_vector = flatten(simulation_eq.time_advance_arrays)
initial_equations += [Eq(conserve_vector[0], d)]
initial_equations += [Eq(conserve_vector[1], d*u0)]
initial_equations += [Eq(conserve_vector[2], d*u1)]
initial_equations += [Eq(conserve_vector[3], d*u2)]
initial_equations += [Eq(conserve_vector[4], p/(gama-1.0) + 0.5* d *(u0**2+u1**2+ u2**2))]
initial_equations += [Eq(GridVariable('temp'), DataObject('x2'))]
initial = GridBasedInitialisation()
initial.add_equations(copy.deepcopy(initial_equations))

multi_block.set_equations([initial])

# block 0 boundary conditions
mb_bcs = {0:None, 1:None, 2:None}
# Boundary conditions for block 0
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and [SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,0] and [2, 1, 0]
block0_bc = []
direction = 0
side = 0
block0_bc.append(InterfaceBC(direction=0, side=0,  match=(1, 0, 0, True)))
block0_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block0_bc.append(SharedInterfaceBC(direction=1, side=0,  match=(2, 1, 0, True)))
block0_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block0_bc.append(PeriodicBC(direction=2, side=0))
block0_bc.append(PeriodicBC(direction=2, side=1))
mb_bcs[0] = block0_bc

# Boundary conditions for block 1
#The boundary conditions are [InterfaceBC, InterfaceBC] in x0 direction and [wall, Inflow]  in x1 direction 
# Matching boundaries are located at are [0,0,0] and [2, 0, 0]
block1_bc = []
block1_bc.append(InterfaceBC(direction=0, side=0,  match=(0, 0, 0, True)))
block1_bc.append(InterfaceBC(direction=0, side=1,  match=(2, 0, 0, False)))
# Wall temperature is required for halo points
Twall = ConstantObject('Twall')
Twall.value = 1.0
wall_energy = [Eq(conserve_vector[-1], Twall*conserve_vector[0]/((gama-1.0)*gama*Minf*Minf))]
block1_bc.append(IsothermalWallBC(direction=1, side=0, equations=wall_energy))
block1_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block1_bc.append(PeriodicBC(direction=2, side=0))
block1_bc.append(PeriodicBC(direction=2, side=1))
mb_bcs[1] = block1_bc

# Boundary conditions for block 2
# The boundary conditions are [InterfaceBC, outflow] in x0 direction and  SharedInterfaceBC, Inflow]  in x1 direction 
# Matching boundaries are located at are [1,0,1] and [0, 1, 0]
block2_bc = []
block2_bc.append(InterfaceBC(direction=0, side=0,  match=(1, 0, 1, False)))
block2_bc.append(ExtrapolationBC(direction=0, side=1, order=0))
block2_bc.append(SharedInterfaceBC(direction=1, side=0,  match=(0, 1, 0, True)))
block2_bc.append(DirichletBC(direction=1, side=1, equations=initial_equations))
block2_bc.append(PeriodicBC(direction=2, side=0))
block2_bc.append(PeriodicBC(direction=2, side=1))
mb_bcs[2] = block2_bc
# Set the multi block boundary conditions
multi_block.set_block_boundaries(mb_bcs)

# Add shock-capturing
# filter_list = []
# for block in multi_block.blocks:
#     ShockFilter = WENOFilter(block, order=3, metrics=metriceq, dissipation_sensor='Constant', Mach_correction=True)
#     filter_list += [ShockFilter.equation_classes]
# multi_block.set_filters(filter_list)

# HDF5 input/output
x,y,z = symbols("x0, x1, x2", **{'cls':DataObject})
kwargs = {'iotype': "Write"}
h5 = iohdf5(save_every=10000, **kwargs)
h5.add_arrays(simulation_eq.time_advance_arrays + [x, y, z])
multi_block.setio([h5])
# Read in the grid file
kwargs = {'iotype': "Read"}
h5_read = iohdf5(**kwargs)
h5_read.add_arrays([x, y, z])
multi_block.setio([h5_read])

# Perform the discretization
multi_block.discretise()
# Add the coordinate excahnges for the multi-block-treatment to the solution of block 2
kernels = create_exchange_calls_codes(multi_block, [x,y])
# Add the wake treatment kernels
wake_ker = generate_wake_kernel(conserve_vector, multi_block, wall_energy[0])
# Sponge kernel for block 0
sponge_ker_block0 = generate_sponge_kernel(conserve_vector, multi_block.get_block(0))
# Sponge kernel for block 2
sponge_ker_block2 = generate_sponge_kernel(conserve_vector, multi_block.get_block(2))

# Add wake exchanges and kernels to block2 boundary conditions
b = multi_block.get_block(2)
for no, eq in enumerate(b.list_of_equation_classes):
    # Add coordinate exchanges to the Block2 GridBasedInitialisation
    if isinstance(eq, GridBasedInitialisation):
        eq.Kernels += kernels
    # Add sponge kernels to block2 spatial solution i.e after evaluating the residuals
    elif isinstance(eq, SimulationEquations):
        eq.Kernels += [sponge_ker_block0, sponge_ker_block2]
        eq.boundary_kernels += wake_ker

# Create the OPS C code
alg = TraditionalAlgorithmRKMB(multi_block)
OPSC(alg)
# NaN check and iteration counter
print_iteration_ops(NaN_check='rho_B0')
# Substitute simulation parameter values
constants = ['gama', 'Minf', 'Pr', 'Re', 'dt', 'niter', 'sigma_filt'] # strength of the free-stream filtering
values = ['1.4', '0.5', '0.72', '50000.0', '0.0001', '100000', '0.01']
# Block 0
constants += ['block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0']
values += ['801', '692', '5', '5.0/(block0np0 - 1.0)', '7.3/(block0np1 - 1.0)', '0.12/block0np2']
# Block 1
constants += ['block1np0', 'block1np1', 'block1np2', 'Delta0block1', 'Delta1block1', 'Delta2block1']
values += ['1799', '692', '5', '7.29705965995/(block1np0 - 1.0)', '7.3/(block1np1 - 1.0)', '0.12/block1np2']
# Block 2
constants += ['block2np0', 'block2np1', 'block2np2', 'Delta0block2', 'Delta1block2', 'Delta2block2']
values += ['801', '692', '5', '5.0/(block2np0 - 1.0)', '7.3/(block2np1 - 1.0)', '0.12/block2np2']
substitute_simulation_parameters(constants, values)

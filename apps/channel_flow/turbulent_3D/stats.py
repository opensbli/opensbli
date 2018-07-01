#!/usr/bin/env python
"""@brief Favre averaged statistics for turbulent channel flow
   @author Satya Pramod Jammy
   @details Current work around for the statistics
"""
from opensbli import *
from opensbli.utilities.user_defined_kernels import *
from opensbli.code_generation.algorithm.common import InTheSimulation, AfterSimulationEnds, BeforeSimulationStarts
# A work around for statistics, later separate class will be provided where both temporal and spatial statistics 

def first_order_moments(ndim, conserve_vector):
    # Create the statistics equations, this shows another way of writing the equations
    # create variables means required
    r_m = symbols('rhomean', **{'cls': DataObject})
    ru_m = symbols('rhou0:%dmean' % ndim, **{'cls': DataObject})

    accumulation_equations = []
    # Accumulate density
    accumulation_equations += [Eq(r_m, r_m + conserve_vector[0])]
    # accumulate_momentum
    for i in range(ndim):
        accumulation_equations += [Eq(ru_m[i], ru_m[i] + conserve_vector[i+1])]

    # Now the normalisation
    normalise_equations = []
    # Dependent on temporal scheme changes FIX
    niter = symbols("niter", **{'cls':ConstantObject})

    normalise_equations += [Eq(r_m, r_m/niter)]
    # accumulate_momentum
    for i in range(ndim):
        normalise_equations += [Eq(ru_m[i], ru_m[i]/niter)]

    return accumulation_equations, normalise_equations

def second_order_moments(ndim, conserve_vector):
    
    # Matrix of favre reynolds stresses i.e rhou_i*u_j if i >= j, only lower triangular matrix is considered
    def RS(i,j):
        if i >= j:
            return DataObject("rhou%du%dmean" % (i, j))
        else:
            return 0
    stresses_reynolds = Matrix(ndim, ndim, RS)
    accumulation = []
    normalisation = []
    # Dependent on temporal scheme changes FIX
    niter = symbols("niter", **{'cls':ConstantObject})
    # creeate a  GridVariable for 1/rho
    rho_inv = GridVariable("rhoinv")
    accumulation += [Eq(rho_inv, 1.0/conserve_vector[0])]
    for i in range(ndim):
        for j in range(ndim):
            if i >= j:
                mean_var = DataObject("rhou%du%dmean" % (i, j))
                comp1 = conserve_vector[i+1]
                comp2 = conserve_vector[j+1]*rho_inv
                accumulation += [Eq(mean_var, comp1 * comp2)]
                normalisation += [Eq(mean_var, mean_var/niter)]
    return accumulation, normalisation

def favre_averaged_stats(ndim, conserve_vector):
    
    accumulation = UserDefinedEquations()
    accumulation.algorithm_place = InTheSimulation()

    normalisation = UserDefinedEquations()
    normalisation.algorithm_place = AfterSimulationEnds()

    first_moms_acc, first_mom_normalise = first_order_moments(ndim, conserve_vector)
    accumulation.add_equations(first_moms_acc)
    normalisation.add_equations(first_mom_normalise)

    sec_moms_acc, sec_mom_normalise = second_order_moments(ndim, conserve_vector)
    accumulation.add_equations(sec_moms_acc)
    normalisation.add_equations(sec_mom_normalise)
    return [accumulation, normalisation]
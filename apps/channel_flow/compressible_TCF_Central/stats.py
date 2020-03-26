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
    niter.datatype = Int()

    normalise_equations += [Eq(r_m, r_m/niter)]
    # accumulate_momentum
    for i in range(ndim):
        normalise_equations += [Eq(ru_m[i], ru_m[i]/niter)]

    E_mean = DataObject('E_mean')
    accumulation_equations += [Eq(E_mean, E_mean + conserve_vector[-1]/conserve_vector[0])]
    normalise_equations += [Eq(E_mean, E_mean/niter)]
    pprint(accumulation_equations)

    return accumulation_equations, normalise_equations

def primitive(ndim, conserve_vector):
    # Create the statistics equations, this shows another way of writing the equations
    # create variables means required
    u_m = symbols('u0:%dmean' % ndim, **{'cls': DataObject})
    niter = symbols("niter", **{'cls':ConstantObject})
    niter.datatype = Int()

    accumulation_equations, normalisation_equations = [], []

    grid_vars = [GridVariable('u%d' % i) for i in range(ndim)]

    a_ins = GridVariable("ains")
    
    p_ins = GridVariable("pins")

    t_ins = GridVariable("tins")

    M_ins = GridVariable("mins")

    mu_ins = GridVariable("muins")

    accumulation_equations += [Eq(grid_vars[i], conserve_vector[i+1]/conserve_vector[0]) for i in range(ndim)]

    pprint(accumulation_equations)
    # accumulate_momentum
    for i in range(ndim):
        accumulation_equations += [Eq(u_m[i], u_m[i] + grid_vars[i])]
    # accumulate_momentum
    for i in range(ndim):
        normalisation_equations += [Eq(u_m[i], u_m[i]/niter)]


    for i in range(ndim):
        for j in range(ndim):
            if i >= j:
                mean_var = DataObject("u%du%dmean" % (i, j))
                comp1 = grid_vars[i]
                comp2 = grid_vars[j]
                accumulation_equations += [Eq(mean_var, (comp1 * comp2) + mean_var)]
                normalisation_equations += [Eq(mean_var, mean_var/niter)]

   
    p_mean = DataObject('p_mean')
    pp_mean = DataObject('pp_mean')
    a_mean = DataObject('a_mean')
    T_mean = DataObject('T_mean')
    TT_mean = DataObject('TT_mean')
    M_mean = DataObject('M_mean')
    mu_mean = DataObject('mu_mean')
    
    p_ins = 0.4 * (conserve_vector[-1] - (0.5*conserve_vector[0]*(grid_vars[0]**2)+0.5*conserve_vector[0]*(grid_vars[1]**2)+0.5*conserve_vector[0]*(grid_vars[2]**2)))

    a_ins = (1.4 * p_ins / conserve_vector[0])**0.5 

    T_ins = p_ins*1.4*0.0955*0.0955/conserve_vector[0]

    mu_ins = T_ins**0.7

    M_ins = ((grid_vars[0]**2)+(grid_vars[1]**2)+(grid_vars[2]**2))**0.5 / a_ins

    accumulation_equations += [Eq(p_mean, p_mean + p_ins)]
    normalisation_equations += [Eq(p_mean, p_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(pp_mean, pp_mean + (p_ins*p_ins))]
    normalisation_equations += [Eq(pp_mean, pp_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(a_mean, a_mean + a_ins)]
    normalisation_equations += [Eq(a_mean, a_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(T_mean, T_mean + T_ins)]
    normalisation_equations += [Eq(T_mean, T_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(TT_mean, TT_mean + (T_ins*T_ins))]
    normalisation_equations += [Eq(TT_mean, TT_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(mu_mean, mu_mean + mu_ins)]
    normalisation_equations += [Eq(mu_mean, mu_mean/niter)]
    pprint(accumulation_equations)

    accumulation_equations += [Eq(M_mean, M_mean + M_ins)]
    normalisation_equations += [Eq(M_mean, M_mean/niter)]
    pprint(accumulation_equations)

    return accumulation_equations, normalisation_equations

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
    niter.datatype = Int()
    
    # creeate a  GridVariable for 1/rho
    rho_inv = GridVariable("rhoinv")
    accumulation += [Eq(rho_inv, 1.0/conserve_vector[0])]
    for i in range(ndim):
        for j in range(ndim):
            if i >= j:
                mean_var = DataObject("rhou%du%dmean" % (i, j))
                comp1 = conserve_vector[i+1]
                comp2 = conserve_vector[j+1]*rho_inv
                accumulation += [Eq(mean_var, (comp1 * comp2) + mean_var)]
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

    prim_acc, prim_norm = primitive(ndim, conserve_vector)
    accumulation.add_equations(prim_acc)
    normalisation.add_equations(prim_norm)

    sec_moms_acc, sec_mom_normalise = second_order_moments(ndim, conserve_vector)
    accumulation.add_equations(sec_moms_acc)
    normalisation.add_equations(sec_mom_normalise)
    return [accumulation, normalisation]

#!/usr/bin/env python

""" Routines for algorithm generation """

from sympy import *
import os

BUILD_DIR = os.getcwd()

import logging
LOG = logging.getLogger(__name__)


class Algorithm(object):

    def __init__(self, temporal_scheme, temporal_order, constant_timestep, spatial_scheme, spatial_order, work_array,
                 evaluation_count, bcs, language, multiblock, explicit_filter):
        self.temporal_scheme = temporal_scheme  # Temporal discretisation scheme (e.g. RK3).
        self.temporal_order = temporal_order  # Temporal order of accuracy
        self.constant_timestep = constant_timestep  # Constant timestep.
        self.spatial_scheme = spatial_scheme  # Spatial discretisation scheme
        self.spatial_order = spatial_order  # Spatial order of accuracy
        self.work_array = work_array + "%d" # Work array
        self.evaluation_count = evaluation_count  # Evaluation count
        self.bcs = bcs  # Boundary conditions
        self.language = language  # Generated code's language
        self.multiblock = multiblock  # Multi-block
        # Multi-block stuff
        if self.multiblock:
            self.nblocks = Symbol('nblocks')
        else:
            self.nblocks = 1
        self.explicit_filter = explicit_filter  # Explicit filter
        
        # Check the input
        self.sanity_check()
        return

    def sanity_check(self):
        """ Check the algorithm options. """

        # If scheme is Runge-Kutta
        if not self.temporal_scheme in ['RK']:
            raise NotImplementedError('Implement %s time stepping scheme' % self.temporal_scheme)

        # Spatial discretisation scheme
        if not self.spatial_scheme in ['central_diff']:
            raise NotImplementedError('Implement %s spatial discretisation scheme' % self.spatial_scheme)

        return

def derivative_formula(system, algorithm):
    ''' finds the finite difference weights for the spatial scheme and order of
	accuracy provided in the algorithm. Updates the fd_weights and fd_points
	in the system

	:args system, algorithm
	:returns none

	'''
    ooa = algorithm.spatial_order
    scheme = algorithm.spatial_scheme
    points = []
    order = 2
    if scheme == 'central_diff':
        points = list(i for i in range(-ooa/2, ooa/2+1))
        if len(points) < order+1:
            raise ValueError("Too few points for derivative of order %d" % order)
        weights = finite_diff_weights(order, points, 0)
        system.fd_weights = weights
        system.fd_points = points
        for dim in range(system.ndim):
            system.halos.append(int(ooa/2))
            system.gridhalo = system.gridhalo + [system.grid[dim+1]+2*system.halos[dim]]
    else:
        raise ValueError("Implement %s scheme in derivative formulation" % scheme)
    return

def sort_evals(inp, evald):
    '''
    Sorts a list of formulas depending on the evaluation i.e. the equations are sorted in
    a way that all the terms in the RHS are evaluated before the equation is evaluated

    :args list of equations and the variables evaluated previously
    :returns evaluated variables and sorted list
    # TODO break if the loop runs for ever
    '''
    inpcopy = [te for te in inp]
    out = []
    while inpcopy:
        for te in inpcopy:
            ter = te.rhs.atoms(Indexed)
            if not ter.difference(set(evald)):
                out.append(te)
                inpcopy.remove(te)
                evald.append(te.lhs)
    return out, evald

def apply_der(der, inp):
    ''' Applies the derivatives formula depending on the finite difference weights
    evaluated in derivative_formula

    : derivative_formula should be called before this
    :args derivative, inputs
    :returns the formula of the derivative in symbolic notation
    '''
    order = len(der.args) - 1
    temp = []
    lwr = inp.block.lower
    upr = inp.block.upper
    for ind in inp.fd_points:
        outderivative = der.args[0]
        for d in der.args[0].atoms(Indexed):
            var = d.base
            index = d.indices
            dire = der.args[1]
            ind1 = str(index[0]).replace('%s' % str(dire), '%s' % str(dire + ind))
            ind1 = Idx(ind1, (lwr, upr))
            outderivative = outderivative.subs(d, var[ind1])
        temp = temp + [outderivative]
    derivative = 0
    N = len(inp.fd_points) - 1
    for nu in range(0, len(inp.fd_points)):
        derivative += inp.fd_weights[order][N][nu]*temp[nu]

    if order == 1:
        derivative = derivative/(Symbol('d%s' % str(dire)))
        inp.constants += [Symbol('d%s' % str(dire))]
    elif order == 2:
        derivative = derivative/(Symbol('d%s' % str(dire))*Symbol('d%s' % str(dire)))
        inp.constants += [Symbol('d%s' % str(dire))]
    else:
        raise NotImplementedError('Implement the derivative of order %d' % order)
    return derivative

#!/usr/bin/env python

""" Routines for algorithm generation """

from sympy import *
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_application)
transformations = standard_transformations + (implicit_application,)
import os

BUILD_DIR = os.getcwd()

import logging
LOG = logging.getLogger(__name__)


class Algorithm(object):

    def __init__(self):
        self.temporal_scheme = []  # Temporal discretisation scheme (e.g. RK3).
        self.temporal_order = []  # Temporal order of accuracy
        self.constant_timestep = True  # Constant timestep.
        self.spatial_scheme = []  # Spatial discretisation scheme
        self.spatial_order = []  # Spatial order of accuracy
        self.work_array = []  # Work array
        self.evaluation_count = []  # Evaluation count
        self.bcs = []  # Boundary conditions
        self.language = []  # Generated code's language
        self.multiblock = False  # Multi-block
        self.nblocks = 1  # Number of blocks
        self.expfilter = []  # Explicit filter
        return

    def read_input(self, path):
        """ Read the algorithm input from the algorithm file.
        
        :arg str path: the path to the algorithm file.
        """
        lines = [line for line in open(path, "r").read().splitlines() if line]

        # Get all the comments in the file
        comment_lineno = []
        for ind, line in enumerate(lines):
            if line[0] == '#':
                comment_lineno.append(ind)

        # Time-stepping scheme
        temp = lines[comment_lineno[0]+1:comment_lineno[1]]
        if len(temp) > 1:
            raise ValueError('Temporal scheme is defined wrongly in algorithm text file')
        temp = temp[0].split(',')
        if temp[0] != 'RK':
            raise ValueError('Implement %s time stepping scheme' % temp[0])
        self.temporal_scheme = temp[0]

        # Temporal order of accuracy
        self.temporal_order = int(temp[1])

        # Constant or local time-step
        temp = lines[comment_lineno[1]+1:comment_lineno[2]]
        if len(temp) > 1:
            raise ValueError('Time step can be const or local only')
        if temp[0] == 'const':
            self.constant_timestep = True
        elif temp[0] == 'local':
            self.constant_timestep = False
        else:
            raise ValueError('Time step can be const or local only')

        # Spatial discretisation scheme
        temp = lines[comment_lineno[2]+1:comment_lineno[3]]
        if len(temp) > 1:
            raise ValueError('Spatial scheme is defined wrongly in algorithm text file')
        temp = temp[0].split(',')
        if temp[0] != 'central_diff':
            raise ValueError('Implement %s spatial scheme' % temp[0])
        self.spatial_scheme = temp[0]

        # Spatial order of accuracy
        self.spatial_order = int(temp[1])

        # Evaluation count
        temp = lines[comment_lineno[3]+1:comment_lineno[4]]
        if len(temp) > 1:
            raise ValueError('The evaluation count provided is wrong')
        self.evaluation_count = int(temp[0])
        temp = lines[comment_lineno[4]+1:comment_lineno[5]]
        if len(temp) > 1:
            raise ValueError('The evaluation count provided is wrong')
        self.work_array = '%s%%d' % temp[0]
        temp = lines[comment_lineno[5]+1:comment_lineno[6]]

        # Read the boundary conditions
        for te in temp:
            te = te.replace(' ', '')
            self.bcs = self.bcs + [tuple(te.strip().split(','))]

        # Language
        temp = lines[comment_lineno[6]+1:comment_lineno[7]]
        self.language = temp[0]

        # Multi-block stuff
        temp = lines[comment_lineno[7]+1:comment_lineno[8]]
        if len(temp) > 1:
            raise ValueError('Blocks defined wrongly in the code')
        temp = temp[0]
        if temp == 'SB':
            self.multiblock = False
            self.nblocks = 1
        elif temp == 'MB':
            self.multiblock = True
            self.nblocks = Symbol('nblocks')
        else:
            raise ValueError('Blocks can be SB for single block or MB for multi-block')

        # Spatial filtering of conservative varibles
        temp = lines[comment_lineno[8]+1:comment_lineno[9]]
        temp = temp[0].split(',')
        for te in temp:
            te = te.strip()
            if te == 'T':
                self.expfilter.append(True)
            elif te == 'F':
                self.expfilter.append(False)
            else:
                raise ValueError('Filter can be T or F only')

        return

def derivative_formula(system, algorithm):
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
    inpcopy = [te for te in inp]
    out = []
    while inpcopy:
        # LOG.debug('in while loop')
        for te in inpcopy:
            ter = te.rhs.atoms(Indexed)
            if not ter.difference(set(evald)):
                out.append(te)
                inpcopy.remove(te)
                evald.append(te.lhs)
    return out, evald

def apply_der(der, inp):
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
        # derivative = derivative.simplify()
        # derivative = ratsimp(derivative)
        # derivative = derivative.simplify()
        # derivative = ratsimp(derivative)
        # outderivative = outderivative.subs(d, derivative)
    if order == 1:
        derivative = derivative/(Symbol('d%s' % str(dire)))
        inp.constants += [Symbol('d%s' % str(dire))]
    elif order == 2:
        derivative = derivative/(Symbol('d%s' % str(dire))*Symbol('d%s' % str(dire)))
        inp.constants += [Symbol('d%s' % str(dire))]
    else:
        raise NotImplementedError('Implement the derivative of order %d' % order)
    # easydebug
    # derivative = derivative.simplify()
    # derivative = ratsimp(derivative)
    # end easydebug
    return derivative

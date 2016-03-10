#!/usr/bin/env python

#    AutoFD: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of AutoFD.

#    AutoFD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    AutoFD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with AutoFD.  If not, see <http://www.gnu.org/licenses/>.

import os
import subprocess
from sympy import *
from sympy.parsing.sympy_parser import *

import logging
LOG = logging.getLogger(__name__)

from .equations import *
from .algorithm import *
from .codegen_utils import *
from .fortran import *
from .opsc import *
from .latex import LatexWriter
from .array import MutableDenseNDimArray
BUILD_DIR = os.getcwd()

class Scheme():
    def __init__(self,scheme, order):
        self.scheme = scheme
        self.order = order
        return

class NumericalGrid():
    def __init__(self, shape):
        self.indices = [Symbol('i%d'%ind, integer = True) for ind, val in enumerate(shape)]
        self.shape = shape
        # as of now uniform grid, later change this as and when required
        self.uniform = [True for ind,val in enumerate(shape)]
        #self.deltas = [EinsteinTerm'delta_i']
        self.time_grid = self.indices + [symbols('Time',integer=True)]
        et = EinsteinTerm('deltai_i');et.is_constant = True
        self.deltas = et.ndimarray(et.IndexedObj(len(shape)))
        self.halos = []
        return
    def work_array(self,name):
        '''
        No shape information will be provided; as the shape of the arrays might change based
        on the computations (including Halos or excluding halos)
        '''
        out = IndexedBase('%s'%name)[self.indices]
        return out

class SpatialDerivative():
    """
    This initializes the spatial derivatives of an arbitrary function 'F'
    on the numerical grid with the provided spatial scheme.
    For wall boundary condition this will have a dependency on grid range need to think of
    that
    """
    def __init__(self, spatial,grid, max_order):
        """
        This initializes the spatial derivative class, which gives the equations
        of spatial Derivatives for combinations of

        Spatial scheme, and order of accuracy
        """
        # stencil should be formula dependant
        self.stencil = [[] for dim in grid.shape]
        self.update_stencil(spatial,grid)
        self.derivatives = []
        self.der_direction = grid.indices
        self.deltas = grid.deltas
        fn = IndexedBase('f',shape = grid.shape)[grid.indices]
        self.fn = fn
        self.Derivative_formulas(fn,max_order, grid)
        return

    def update_stencil(self,spatial, grid):
        for dim, val in enumerate(grid.shape):
            if spatial.scheme == 'central':
                points = list(i for i in range(-spatial.order/2, spatial.order/2+1))
                grid.halos.append(tuple([-spatial.order/2, spatial.order/2]))
            else:
                raise NotImplementedError("Only central difference schemes are supported")
            self.stencil[dim] = [grid.indices[dim] + i for i in points]
            #print points, grid.indices[dim], self.stencil[dim]
        return
    def Derivative_formulas(self, fn, max_order, grid):
        derivatives = []
        derivatives += [fn] # later change this to interpolation
        derivative_formula = []
        derivative_formula += [fn]
        comp_kernels = []
        comp_kernels += [fn]
        for order in range(1,max_order+1):
            shape = tuple( [len(grid.indices) for ind in range(order)])
            array = MutableDenseNDimArray.zeros(*shape)
            fdarray = MutableDenseNDimArray.zeros(*shape)
            deriv_kernel = MutableDenseNDimArray.zeros(*shape)
            for ind in np.ndindex(*array.shape):
                der_args = [grid.indices[i] for i in ind]
                name = [str(arg) for arg in ind]
                #name = tuple(ind)
                name = "[%d][%s]"%(order,','.join(name))
                deriv_kernel[ind] = Symbol(name)
                array[ind] = fn.diff(*der_args)
                # find the finite difference formula
                if order == 1 or len(set(der_args)) ==1:
                    fdarray[ind] = as_finite_diff(array[ind], self.stencil[ind[0]])*pow(grid.deltas[ind[0]],-order)
                else:
                    newder = array[ind].subs(derivatives[order-1][ind[:-1]],derivative_formula[order-1][ind[:-1]])
                    fdarray[ind] = as_finite_diff(newder,self.stencil[ind[-1]], wrt=grid.indices[ind[-1]])*pow(grid.deltas[ind[-1]],-1)
            derivatives.append(array)
            derivative_formula.append(fdarray)
            comp_kernels.append(deriv_kernel)
        self.derivatives = derivatives
        self.derivative_formula = derivative_formula
        self.deriv_kernel = comp_kernels
        return
    def get_derivativeformula(self, derivative, order):
        '''This returns the formula for the derivative using the functions provided
        for getting a symbolic derivative for a general function use get_derivative
        used for ceval stuff
        '''
        order = len(derivative.args[1:])
        inds = []
        for arg in derivative.args[1:]:
            inds = inds + [self.der_direction.index(arg)]
        if order == 1 or len(set(inds)) ==1:
            formula = as_finite_diff(derivative, self.stencil[inds[0]])*pow(self.deltas[inds[0]],-order)
        else:
            loweder = Derivative(derivative.args[0], *inds[:-1])
            raise ValueError("first update the derivative of %s before calling %s"%(loweder, derivative))
        return formula
    def get_derivative(self, derivative):
        '''
        This returns a tuple to which the derivaitve formula exists in
        already evaluated derivatives.
        '''
        order = len(derivative.args[1:])
        inds = []
        for arg in derivative.args[1:]:
            inds = inds + [self.der_direction.index(arg)]
        generalformula = []
        subevals = []
        requires = []
        if order == 1 or len(set(inds)) ==1:
            generalformula += [order,tuple(inds)]
            if len(derivative.args[0].atoms(Indexed)) >1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += [derivative.args[0]]
        else:
            if len(derivative.args[0].atoms(Indexed)) >1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))

            else:
                subevals += [None]
                requires += [derivative.args[0]]
            generalformula += [order-1,tuple([inds[-1]])]
            requires += [self.derivatives[order-1][inds[:-1]].subs(self.fn,derivative.args[0])]

        return generalformula, subevals, requires


class TimeDerivative():
    def __init__(self, temporal,grid, const_dt):
        # This can be changed totally but as of now leave it
        if const_dt:
            dt = EinsteinTerm('deltat')
            dt.is_constant = True; dt.is_commutative = True
        else:
            raise NotImplementedError("Varying delta t is not implemented in the code")
        fn = IndexedBase('f')[grid.time_grid]
        time_index = grid.time_grid[-1]
        self.coeffs = {}
        self.timeder = fn.diff(time_index)
        self.formula = []
        self.stages = 0
        if temporal.scheme == "Forward":
            if temporal.order ==1:
                self.formula = as_finite_diff(self.timeder, [time_index+1, time_index])*pow(dt,-1)
                self.stages = 1
                return
            else:
                raise ValueError("Forward time stepping scheme of first order accuracy is allowed")
        elif temporal.scheme == "RungeKutta":
            raise NotImplementedError("IMPLEMENT RUNGE KUTTA TIME STEPPING SCHEME")
        else:
            raise ValueError("Only forward or Runge-Kutta time stepping schemes are allowed")
        return
def maximum_derivative_order(equations):
    order = set()
    for eq in equations:
        for atom in eq.atoms(Derivative):
            order.add(len(atom.args)-1)
    return max(order)
def vartoGridArray(variable,grid):
    '''
    Converts a variable/ function or Indexed Object to a Indexed base on the Grid
    inputs: variable, grid
    returns: the Grid array
    '''
    if isinstance(variable, Indexed):
        return variable.base[grid.indices]
    elif isinstance(variable, Function):
        return IndexedBase('%s'%variable.func)[grid.indices]
    else:
        raise ValueError("Only functions or Indexed Objects are supported", variable)
    return
def get_spatial_derivatives(equations):
    ders = []
    count = {}
    time_ders = []
    for eq in equations:
        pot = preorder_traversal(eq)

        for p in pot:
            if p in ders:
                pot.skip()
                count[p] = count[p]+1
                continue
            elif isinstance(p, Derivative):
                if all(arg != EinsteinTerm('t') for arg in p.args):
                    pot.skip()
                    ders.append(p)
                    count[p] = 1
                else:
                    pot.skip()
                    time_ders.append(p)
            else:
                continue
    return ders, count, time_ders
def get_grid_variables(equations):
    variables = []
    count = {}
    for eq in equations:
        pot = preorder_traversal(eq)
        for p in pot:
            if p in variables:
                pot.skip()
                count[p] = count[p]+1
                continue
            elif isinstance(p, Indexed):
                pot.skip()
                variables.append(p)
                count[p] = 1
            else:
                continue
    return variables, count

class Evaluations():
    def __init__(self, lhs,rhs, requires,subevals = None, wk=None):
        self.store = True
        if isinstance(lhs, Derivative):
            self.is_der = True
            self.is_formula = False
            if subevals:
                self.subevals = subevals
            else:
                self.subevals = [None]
            if wk:
                self.work = wk
            else:
                self.work = None
            self.formula = rhs
            self.requires = requires
            self.evaluation_range = []
        else:
            self.is_formula = True
            self.is_der = False
            if subevals:
                self.subevals = subevals
            else:
                self.subevals = [None]
            if wk:
                self.work = wk
            else:
                self.work = None
            self.formula = rhs
            self.requires = requires
            self.evaluation_range = []
        return
class ComputationalKernel():
    def __init__(self, equations, ranges, SD=None):
        '''
        This should do two things
        1. Able to write the kernel calling function based on language
            For this we require (IN OPSC)
            a. Name of the kernel to call
            b. Block on which it should work, dimensions of the block
            c. ins, outs, inouts and their stencil of access, data type
            d. Indices of the array if required
        2. Write the computational kernel, requires writing kernel header
        and the computations
            Require for OPSC
            a.Kernel header can be written with ins, outs,
        inouts and indices along with the data type
            b. TO write the computations the Indexed Objects are to be modified
        this modification requires OPS_ACC values which can also be populated and
        grid indices are replaced with 0's
        All in all we require
        1. Name
        2. Block (updated from the call to ops_write)
        3. ins, outs, inouts and so on
        4. Stencils of access
        '''
        self.ranges = ranges# range of the kernel
        self.name = None # None generates automatic kernel name
        if isinstance(equations, list):
            self.equations = equations
        else:
            self.equations = [equations]

        self.inputs = {}
        self.outputs = {}
        self.inputoutoutput = {}
        self.constants = {}
        if not SD:
            self.classify_indexed_objects()
        else:
            pass
        return
    def classify_indexed_objects(self):
        ins = []
        outs = []
        inouts = []
        consts = []
        for eq in self.equations:
            ins = ins + list(eq.rhs.atoms(Indexed))
            outs = outs + list(eq.lhs.atoms(Indexed))
            consts = consts + list(eq.atoms(EinsteinTerm))
        inouts = set(outs).intersection(set(ins))
        ins = set(ins).difference(inouts)
        outs = set(outs).difference(inouts)
        indexbase_ins = set([v.base for v in ins])
        indexbase_outs = set([v.base for v in outs])
        indexbase_inouts = set([v.base for v in inouts])
        for v in indexbase_ins:
            indexes = [vin.indices for vin in ins if vin.base==v]
            self.inputs[v] = indexes
        for v in indexbase_outs:
            indexes = [vin.indices for vin in outs if vin.base==v]
            self.outputs[v] = indexes
        for v in indexbase_inouts:
            indexes = [vin.indices for vin in inouts if vin.base==v]
            self.inputoutoutput[v] = indexes
        # TODO Add any Idx objects to the inputs as well if have one
        consts = set(consts)
        return
    def classify_ders(self):

        return
class SpatialSolution():
    def __init__(self, equations, formulas, grid, spatial_scheme):
        alleqs = flatten(list(e.expanded for e in equations))
        allformulas = flatten(list(e.expanded for e in formulas))
        max_order = maximum_derivative_order(alleqs)
        SD = SpatialDerivative(spatial_scheme,grid,max_order)
        grid_arrays = {}
        range_used = {}
        grid_variables, variable_count = get_grid_variables(alleqs+allformulas)
        for atom in grid_variables:
            grid_arrays[atom] = vartoGridArray(atom, grid)
        spatialders, dercount, time_derivatives = get_spatial_derivatives(alleqs+allformulas)
        # Define the formulas on the grid, this is substituting the old with new
        # TODO do a sanity check of the formulas, i.e. remove all the formulas that
        # are not used in the equations
        evals = {}
        for form in allformulas:
            out = form
            for atom in form.atoms(Indexed):
                out = out.subs(atom, grid_arrays[atom])
            evaluated = Evaluations(out.lhs,out.rhs, list(out.rhs.atoms(Indexed)), None,out.lhs)
            evals[out.lhs] = evaluated
        et = EinsteinTerm('x_i');et.is_constant = True
        coord = et.ndimarray(et.IndexedObj(len(grid.shape)))
        coord = coord.tolist()
        wkarray = 'wk'; wkind = 0;
        for der in spatialders:
            out = der# Modify the derivative to be a derivative on grid
            wk = grid.work_array('%s%d'%(wkarray,wkind)); wkind = wkind +1
            for atom in der.atoms(Indexed):
                out = out.subs(atom, grid_arrays[atom])
            for arg in out.args[1:]:
                out = out.subs(arg, grid.indices[coord.index(arg)])
            generalformula, subevals, requires = SD.get_derivative(out)
            grid_arrays[der] = out
            evaluated = Evaluations(out,generalformula,requires, subevals, wk)
            evals[out] = evaluated
        # we will assume that all the functions in time derivative are known at the start
        known = [grid_arrays[d.args[0]] for d in time_derivatives]
        for val in known:
            evaluated = Evaluations(val,val, None, None,val)
            evals[val] =  evaluated
        # Sort the Formulas
        orderof_evaluations = [grid_arrays[d.args[0]] for d in time_derivatives]
        orderof_evaluations = sort_evaluations(orderof_evaluations,evals, Indexed)
        # sort the derivatives
        orderof_evaluations = sort_evaluations(orderof_evaluations,evals, Derivative)
        # update the range of evaluations for each evaluation
        range_of_evaluation(orderof_evaluations, evals,grid, SD)
        # now define a ComputationalKernel for each of the evaluations

        ''' All the primitive variables (IndexedObjects)
        excluding those which have a time derivative are stored into a kernel
        '''
        forms = [ev for ev in orderof_evaluations if isinstance(ev, Indexed) and ev not in known]
        ranges = [evals[ev].evaluation_range for ev in forms]
        subevals = flatten([evals[ev].subevals for ev in forms])
        subeval_truth = [ev == None for ev in subevals]
        # check if all the ranges of evaluation are the same for the Formulas
        range_truth = [ranges[0][i] == val[i] for val in ranges for i in range(len(ranges[0]))]
        computations = []
        eqs = []
        eqs = [Eq(evals[ev].work, evals[ev].formula) for ev in forms]

        # if same range then combine them into a single computation else store into different computations
        if all(range_truth) and all(subeval_truth):
            computations.append(ComputationalKernel(eqs, ranges))
        else:
            for number,eq in enumerate(eqs):
                computations.append(ComputationalKernel(eq, ranges[number]))
        # Now process the Derivatives
        derivatives = [ev for ev in orderof_evaluations if isinstance(ev, Derivative) and ev not in known]
        ranges = [evals[ev].evaluation_range for ev in derivatives]
        subevals = [evals[ev].subevals for ev in derivatives]
        require = [evals[ev].requires for ev in derivatives]
        for number,der in enumerate(derivatives):
            if not any(isinstance(req, Derivative) for req in require[number]):
                if all(subev == None for subev in subevals[number]):
                    rhs = SD.get_derivativeformula(der,evals[der].formula[0])
                    eq = Eq(evals[der].work,rhs)
                    computations.append(ComputationalKernel(eq, ranges[number]))
                    #print "Computed"
                else:
                    # store into temporary array sub evaluation
                    eqs = []
                    tempwkind = wkind
                    for subev in subevals[number]:
                        wk = grid.work_array('%s%d'%(wkarray,tempwkind));tempwkind = tempwkind +1
                        for req in require[number]:
                            local_range = [req,evals[req].evaluation_range]
                            subev = subev.subs(req,evals[req].work)
                        eqs.append(Eq(wk,subev))
                    computations.append(ComputationalKernel(eqs, local_range))
                    for eq in eqs:
                        newder = der.subs(eq.rhs,eq.lhs)
                    rhs = SD.get_derivativeformula(newder,evals[der].formula[0])
                    eq = Eq(evals[der].work,rhs)
                    computations.append(ComputationalKernel(eq, ranges[number]))
            else:
                newder = der
                if all(subev == None for subev in subevals[number]):
                    for req in require[number]:
                        newder = newder.subs(req,evals[req].work)
                else:
                    raise NotImplementedError("Sub evaluations in a mixed derivative")
                rhs = SD.get_derivativeformula(newder,evals[der].formula[0])
                eq = Eq(evals[der].work,rhs)
                computations.append(ComputationalKernel(eq, ranges[number]))
        updated_eq = [eq for eq in alleqs]
        # All the spatial computations are performed now get the updated equations
        for eqno,eq in enumerate(updated_eq):
            spatialders, dercount, time_derivatives = get_spatial_derivatives([eq])
            grid_variables, variable_count = get_grid_variables([eq])
            spatialders = (sorted(spatialders, cmp = decreasing_order))
            # substitute spatial ders first
            for var in spatialders+grid_variables:
                new = evals[grid_arrays[var]].work
                updated_eq[eqno] = updated_eq[eqno].subs(var,new)
            pprint(updated_eq[eqno])
        return

def range_of_evaluation(orderof_evaluations, evaluations, grid, sdclass):
    '''
    First the ranges of derivatives are updated, then other ranges are updated
    '''
    ders = []
    for ev in orderof_evaluations:
        if isinstance(ev, Derivative):
            ders.append(ev)
        evaluations[ev].evaluation_range = [tuple([0,s]) for s in grid.shape]
    # Update the range for the derivatives
    grouped_ders = group_derivatives(ders)
    for key, value in grouped_ders.iteritems():
        for val in value:
            require = evaluations[val].requires
            form  = evaluations[val].formula
            dire = form[1][0]
            halos = grid.halos[dire]
            for req in require:
                erange = list(evaluations[req].evaluation_range[dire])
                if erange[0] == 0 and erange[1] == grid.shape[dire]:
                    erange[0] = erange[0]+halos[0]
                    erange[1] = erange[1]+halos[1]
                evaluations[req].evaluation_range[dire] = tuple(erange)
    #update the range for the formulas
    for ev in orderof_evaluations:
        if isinstance(ev, Indexed):
            require = evaluations[ev].requires
            if require:
                for req in require:
                    evaluations[req].evaluation_range = evaluations[ev].evaluation_range

    return

def sort_evaluations(evaluated, evaluations, typef):
    for key in evaluations.keys():
        if isinstance(key, typef) and not key in evaluated:
            if all(ev in evaluated for ev in evaluations[key].requires):
                evaluated.append(key)
            else:
                for val in evaluations[key].requires:
                    if not val in evaluated:
                        sort_evaluations(evaluated, {val:evaluations[val]}, typef)
                evaluated.append(key)
    return evaluated
def decreasing_order(s1, s2):
    return cmp(len(s2.args), len(s1.args))
def increasing_order(s1, s2):
    return cmp(len(s1.args), len(s2.args))
def group_derivatives(spatialders):
    spatial_der_dict ={}
    for der in spatialders:
        if der.args[0] in spatial_der_dict.keys():
            spatial_der_dict[der.args[0]] += [der]
        else:
            spatial_der_dict[der.args[0]] = [der]

    for key,value in spatial_der_dict.iteritems():
        if len(value)>1:
            spatial_der_dict[key] = (sorted(value, cmp=increasing_order))
    return spatial_der_dict
class System(object):
    '''
    Re arrange system with new spatial derivatives class etc..
    '''

    def __init__(self):
        return

    def prepare(self, equations, formulas, algorithm, simulation_parameters, initial_conditions):
        expanded_equations = flatten(list(e.expanded for e in equations))
        expanded_formulas = flatten(list(f.expanded for f in formulas))
        variables = flatten(list(e.variables for e in equations))
        conser = flatten(list(e.conser for e in equations))
        constants = list(set(flatten(list(e.constants for e in equations))))
        all_variables = list(set(conser + variables))

        self.ndim = equations[0].ndim
        work_array_index = 0
        self.constants = []
        self.dats = []
        self.grid = []
        # Race_check
        self.race_check = True

        if algorithm.language == 'OPSC':
            l = 0
        elif algorithm.language == 'F90':
            l = 1

        indices = []
        for i in range(self.ndim):
            indices += ['x%d' % i]
        indices = ','.join(indices)
        variable_indices = Idx('%s' % (indices), (l, Symbol('nblock', integer=True)))

        # Single-block and multi-block
        self.block = Idx('blk', (l, Symbol('nblock', integer=True)))
        self.grid = self.grid + [self.block.upper]
        self.nblocks = algorithm.nblocks
        self.multiblock = algorithm.multiblock
        self.blockdims = []
        self.halos = []
        self.gridhalo = []
        self.ranges = {}
        self.iterrange = 0
        self.kernel_index = 0

        if algorithm.language == 'OPSC':
            if self.multiblock:
                raise NotImplementedError("Implement Multi-Block code for OPSC")
            code_file = open(BUILD_DIR+'/%s.cpp' % simulation_parameters["name"], 'w')
            kernel_file = open(BUILD_DIR+'/auto_kernel.h', 'w')
            self.stencils = {}
            self.stencil_index = 0
            self.stencil_name = 'stencil%dD_%%d' % self.ndim
            self.block_name = 'auto_block_OPSC'
            self.kernel_name = 'auto_kernel_%d'
            for dim in range(self.ndim):
                temp = Idx('i%d' % dim, Symbol('nx%dp[blk]' % dim, integer=True))
                self.blockdims.append(temp)
                self.grid = self.grid + [temp.upper+1]
        elif algorithm.language == 'F90':
            if self.multiblock:
                raise NotImplementedError("Implement Multi-Block code for F90")
            code_file = open(BUILD_DIR+'/%s.f90' % simulation_parameters["name"], 'w')
            kernel_file = open(BUILD_DIR+'/subroutines.f90', 'w')
            module_file = open(BUILD_DIR+'/param_mod.f90', 'w')
            self.module = []
            for dim in range(self.ndim):
                temp = Idx('i%d' % dim, (1, Symbol('nx%dp' % dim, integer=True)))
                self.blockdims.append(temp)
                self.grid = self.grid + [temp.upper]
            self.kernel_name = 'subroutine_%d'
        else:
            raise NotImplementedError('Implement indexing for language %s' % (algorithm.language))

        # Check that all the system's properties are set up according to the algorithm.
        self.sanity_check(algorithm)

        # Create new indexed variables and substitute these into the equations in place of the (non-indexed) variables.
        indexed_variables = []
        for variable in all_variables:
            new = Indexed('%s' % str(variable), variable_indices)
            indexed_variables = indexed_variables + [Indexed('%s' % str(variable), variable_indices)]
            for e in range(len(expanded_equations)):
                expanded_equations[e] = expanded_equations[e].subs(variable, new)

        # Prepare the formulas by indexing their variables.
        formula_variables = flatten(list(f.conser for f in formulas))
        formula_variables += flatten(list(f.variables for f in formulas))
        formula_variables = list(set(formula_variables))
        self.constants += list(set(flatten(list(f.constants for f in formulas))))
        for variable in formula_variables:
            new = Indexed('%s' % str(variable), variable_indices)
            for f in range(len(expanded_formulas)):
                expanded_formulas[f] = expanded_formulas[f].subs(variable, new)
        expanded_formulas_dict = equations_to_dict(expanded_formulas)

        # pprint(expanded_equations)
        indexed_variables = set(indexed_variables)
        substis = {}
        evaluated = []
        sorteval = []
        self.conser = []
        for con in conser:
            new = Indexed('%s' % str(con), variable_indices)
            evaluated = evaluated + [new]
            sorteval = sorteval + [new]
            self.conser = self.conser + [new]
            self.dats = self.dats + [new]
        indexed_variables = indexed_variables.difference(set(evaluated))

        # Finished preparing all the stuff now prepare the equations, now the evaluations of primitive variables
        formula_evaluations = []

        for v in indexed_variables:
            val = expanded_formulas_dict.get(v)
            count = variable_count(v, expanded_equations)
            count = algorithm.evaluation_count
            if count >= algorithm.evaluation_count:
                if val:
                    formula_evaluations.append(Eq(v, val))
                    self.dats = self.dats + [v]
                else:
                    raise ValueError('I dont know the formula for %s ' % v)
            else:
                substis[v] = val
                # raise ValueError('Implement how to do for count > evaluation_count')

        # Runge-Kutta timestepping scheme
        if algorithm.temporal_scheme == 'RK':
            rkloop = Idx('nrk', (l, algorithm.temporal_order))
            time_loop = Idx('iter', (l, Symbol('niter', integer=True)))
            a1 = Indexed('a1', Idx('nrk', (l, algorithm.temporal_order)))
            a2 = Indexed('a2', Idx('nrk', (l, algorithm.temporal_order)))
            rk = Indexed('rk', Idx('nrk', (l, algorithm.temporal_order)))
            self.constants += [a1, a2, time_loop.upper+1]
            saveeq = []
            conser_old = []

        formula_evaluations, evaluated = sort_evals(formula_evaluations, evaluated)
        derivatives = flatten(list(e.atoms(Derivative) for e in expanded_equations))
        derivatives = list(set(derivatives))
        derivative_evaluations = []

        # Get the coefficients of the finite difference scheme
        derivative_formula(self, algorithm)
        tempder = []
        for no, der in enumerate(derivatives):
            tempder = tempder + [der]
            if der.args[0] in evaluated:
                val = der.args[0]
            else:
                val = substis.get(der.args[0])
            if algorithm.evaluation_count == 0:
                count = 100
            elif algorithm.evaluation_count == 1000:
                count = 999
            else:
                count = variable_count(der, expanded_equations)
            order = len(der.args) - 1
            if val:
                if order == 1 and count >= algorithm.evaluation_count:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(var, self)
                            derivative_evaluations.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                    else:
                        if algorithm.temporal_scheme == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), variable_indices)
                            if algorithm.constant_timestep:
                                time = Symbol('dt')
                                self.constants += [time]
                            else:
                                time = Indexed('dt', variable_indices)
                                self.dats = self.dats + [time]
                            rkfor = (con - old)/(rk*time)
                            conser_old.append(old)
                            self.dats = self.dats + [old]
                            substis[der] = rkfor
                            saveeq.append(Eq(old, con))

                    # substitute RK3 routine now
                elif order == 2 and count >= algorithm.evaluation_count:
                    if der.args[1] == der.args[2]:
                        var = der.subs(der.args[0], val)
                        temp = algorithm.work_array % (work_array_index)
                        work_array_index += 1
                        new = Indexed('%s' % str(temp), variable_indices)
                        self.dats = self.dats + [new]
                        derf = apply_der(var, self)
                        derivative_evaluations.append(Eq(new, derf))
                        substis[der] = new
                    else:
                        expr = Derivative(der.args[0], der.args[2])
                        if substis.get(expr):
                            var = der.subs(expr, substis.get(expr))
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(var, self)
                            derivative_evaluations.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                        else:
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(expr, self)
                            substis[expr] = new
                            derivative_evaluations.append(Eq(new, derf))
                            self.dats = self.dats + [new]
                            # substis[expr] = new
                            var = der.subs(expr, new)
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(var, self)
                            derivative_evaluations.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                elif order == 1 and count < algorithm.evaluation_count:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            derf = apply_der(var, self)
                            substis[der] = derf
                    else:
                        if algorithm.temporal_scheme == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), variable_indices)
                            if algorithm.constant_timestep:
                                time = Symbol('dt')
                                self.constants += [time]
                            else:
                                time = Indexed('dt', variable_indices)
                                self.dats = self.dats + [time]
                            rkfor = (con - old)/(rk*time)
                            conser_old.append(old)
                            self.dats = self.dats + [old]
                            substis[der] = rkfor
                            saveeq.append(Eq(old, con))
                elif order == 2 and count < algorithm.evaluation_count:
                    if der.args[1] == der.args[2]:
                        var = der.subs(der.args[0], val)
                        derf = apply_der(var, self)
                        substis[der] = derf
                    else:
                        expr = Derivative(der.args[0], der.args[2])
                        if substis.get(expr):
                            var = der.subs(expr, substis.get(expr))
                            derf = apply_der(var, self)
                            # easy_debug
                            # derf = derf.simplify()
                            # easy_debug
                            substis[der] = derf
                        else:
                            derf = apply_der(expr, self)
                            var = der.subs(expr, derf)
                            derf = apply_der(var, self)
                            # for atom in list(derf.atoms(Indexed)):
                                # var = der.subs(expr,atom)
                                # derivative = apply_der(var,self)
                                # derf = derf.subs(atom,derivative)
                            # easy_debug
                            # derf = derf.simplify()
                            # easy_debug
                            substis[der] = derf
                else:
                    raise ValueError('No know;')

            elif count >= algorithm.evaluation_count:
                temp = algorithm.work_array % (work_array_index)
                work_array_index += 1
                new = Indexed('%s' % str(temp), variable_indices)
                self.dats = self.dats + [new]
                var = der.args[0]
                formula_evaluations.append(Eq(new, var))
                # substis[var] = new
                var = der.subs(var, new)
                temp = algorithm.work_array % (work_array_index)
                work_array_index += 1
                new = Indexed('%s' % str(temp), variable_indices)
                derf = apply_der(var, self)
                derivative_evaluations.append(Eq(new, derf))
                substis[der] = new
                self.dats = self.dats + [new]
            elif count < algorithm.evaluation_count:
                derf = apply_der(der, self)
                substis[der] = derf
            else:
                pprint(var)
                raise ValueError('This is wrong')

        evaluations = {}
        evaluation_index = 0
        formula_evaluations, sorteval = sort_evals(formula_evaluations, sorteval)
        evaluations[evaluation_index] = formula_evaluations
        evaluation_index += 1
        for der in derivative_evaluations:
            evaluations[evaluation_index] = der
            evaluation_index += 1

        # Write out algorithm in LaTeX form
        latex = LatexWriter()
        latex.open(path=BUILD_DIR+'/algorithm.tex')
        metadata = {"title": "Algorithm", "author": "Satya P Jammy", "institution": "University of Southampton"}
        latex.write_header(metadata)

        latex.write_string("The equations are\n\n")
        latex.write_equations(expanded_equations, variable_indices)
        latex.write_string('The save state equations are\n\n')
        latex.write_equations(saveeq, variable_indices)
        latex.write_string('The evaluations performed are\n\n')
        latex.write_equations(evaluations, variable_indices)
        latex.write_string('The substitutions in the equations are\n\n')
        latex.write_equations(substis, variable_indices)

        for key, value in substis.iteritems():
            for eqno in range(len(expanded_equations)):
                expanded_equations[eqno] = expanded_equations[eqno].xreplace({key: value})

        # Get the final Runge-Kutta update equations
        tempdict = equations_to_dict(saveeq)
        savedic = dict(zip(tempdict.values(), tempdict.keys()))
        # for race_check errors
        # 1. Find the residue,
        # 2. Update conserve
        # 3. Update old conservative so that there will be no race_errors
        if self.race_check:
            final_equations = {}
            resdue_eq = []
            upd_conser = []
            upd_old_cons = []
            race_ind = 0
        else:
            final_equations = []
        for eqno in range(len(expanded_equations)):
            lh = expanded_equations[eqno].lhs
            rh = expanded_equations[eqno].rhs
            if self.race_check:
                race_ind = race_ind + 1
                temp = 'race_eq%d' % race_ind
                new = Indexed('%s' % str(temp), variable_indices)
                self.dats = self.dats + [new]
                resdue_eq = resdue_eq + [Eq(new, rh)]
                temp = solve(Eq(lh, 0), self.conser[eqno])
                nrdr = fraction(lh)
                eqn = nrdr[1]*new + temp[0]
                eqn1 = eqn.xreplace({rk: a1})
                upd_conser = upd_conser + [Eq(self.conser[eqno], eqn1)]
                old = savedic.get(self.conser[eqno])
                eqn1 = eqn.xreplace({rk: a2})
                upd_old_cons = upd_old_cons + [Eq(old, eqn1)]
            else:
                eqn = Symbol('eqn%d' % eqno)
                final_equations = final_equations + [Eq(eqn, rh)]
                nrdr = fraction(lh)
                temp = solve(Eq(lh, 0), self.conser[eqno])
                eqn = nrdr[1]*eqn + temp[0]
                eqn1 = eqn.xreplace({rk: a1})
                final_equations = final_equations + [Eq(self.conser[eqno], eqn1)]
                old = savedic.get(self.conser[eqno])
                eqn1 = eqn.xreplace({rk: a2})
                final_equations = final_equations + [Eq(old, eqn1)]

        # Apply filtering to the equations
        if self.race_check:
            race_ind = 0
            final_equations[race_ind] = resdue_eq
            race_ind = race_ind+1
            final_equations[race_ind] = upd_conser
            race_ind = race_ind+1
            final_equations[race_ind] = upd_old_cons
            race_ind = race_ind+1

        latex.write_string('The equations after substitutions are\n\n')
        latex.write_equations(expanded_equations, variable_indices)

        latex.write_string('The final rk3 update equations are\n\n')
        latex.write_equations(final_equations, variable_indices)

        if any(algorithm.explicit_filter):
            self.filtereqs = expfiltering(algorithm, self, savedic)
            latex.write_string('The filter equations are\n\n')
            for e in self.filtereqs:
                latex.write_equations(e, variable_indices)

        algorithm_template = {'a00': 'header', 'a01': 'defdec', 'a02': 'grid', 'a03': 'init', 'a04': 'Bcst'}
        algorithm_template['a05'] = 'timeloop'
        algorithm_template['a06'] = 'saveeq'
        algorithm_template['a07'] = 'innerloop'
        algorithm_template['a08'] = 'solveeq'
        algorithm_template['a09'] = 'Bcs'
        algorithm_template['a10'] = 'endinnerloop'
        algorithm_template['a11'] = 'afterinnerloop'
        algorithm_template['a12'] = 'endtimeloop'
        algorithm_template['a13'] = 'aftertimeloop'
        algorithm_template['a14'] = 'footer'

        final_algorithm = {}
        self.constants = list(set(self.constants))
        self.dats = list(set(self.dats))

        # Writing kernels and kernel calls
        if algorithm.language == 'OPSC':
            opsc = OPSC()
            calls = []
            kernels = []
            if evaluations:
                call, kernel = opsc.generate(evaluations, self)  # Evaluations
                calls = calls + call
                kernels = kernels + kernel

            call, kernel = opsc.generate(final_equations, self)  # Evaluations
            calls = calls + call
            kernels = kernels + kernel
            final_algorithm[algorithm_template.get('a08')] = calls

            call, kernel = opsc.generate(saveeq, self)
            kernels = kernels + kernel
            final_algorithm[algorithm_template.get('a06')] = call

            # Initial conditions kernel
            latex.write_equations(initial_conditions, variable_indices)
            call, kernel = opsc.generate(initial_conditions, self)
            kernels = kernels + kernel
            final_algorithm[algorithm_template.get('a03')] = call

        elif algorithm.language == 'F90':
            calls = []
            kernels = []
            if evaluations:
                call, kernel = F90_write_kernel(evaluations, self, algorithm)  # Evaluations
                calls = calls + call
                kernels = kernels + kernel
            call, kernel = F90_write_kernel(final_equations, self, algorithm)  # Evaluations
            calls = calls + call
            kernels = kernels + kernel

            final_algorithm[algorithm_template.get('a08')] = calls

            call, kernel = F90_write_kernel(saveeq, self, algorithm)
            kernels = kernels + kernel
            final_algorithm[algorithm_template.get('a06')] = call
            init_eq = []
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], parse_expr('sin(x)* cos(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], parse_expr('-cos(x)* sin(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            latex.write_equations(init_eq, variable_indices)
            call, kernel = F90_write_kernel(init_eq, self, algorithm)
            kernels = kernels + kernel
            final_algorithm[algorithm_template.get('a03')] = call

        # TODO: Apply boundary conditions
        if algorithm.language == 'OPSC':
            bcs(self, algorithm)
        elif algorithm.language == 'F90':
            self.bccall = ['']

        final_algorithm[algorithm_template.get('a00')] = header_code(self, algorithm.language)
        final_algorithm[algorithm_template.get('a01')] = defdec(self, algorithm, simulation_parameters)
        final_algorithm[algorithm_template['a04']] = self.bccall
        final_algorithm[algorithm_template['a09']] = self.bccall

        temp = '%s Write the grid here if required' % COMMENT_DELIMITER[algorithm.language]
        temp = temp + '\n\n\n' + '%s Grid writing ends here' % COMMENT_DELIMITER[algorithm.language]
        final_algorithm[algorithm_template.get('a02')] = [temp]
        final_algorithm[algorithm_template.get('a11')] = ['%s after  after inner time loop' % COMMENT_DELIMITER[algorithm.language]]

        tloop, tenloop = loop([time_loop], algorithm)
        final_algorithm[algorithm_template.get('a05')] = tloop
        final_algorithm[algorithm_template.get('a12')] = tenloop
        tloop, tenloop = loop([rkloop], algorithm)

        final_algorithm[algorithm_template.get('a07')] = tloop
        final_algorithm[algorithm_template.get('a10')] = tenloop

        # Write the conservative variables to files (including halos)
        final_algorithm[algorithm_template.get('a13')] = after_time(self, algorithm.language)
        final_algorithm[algorithm_template.get('a14')] = footer_code(self, algorithm.language)

        final_algorithm['kernels'] = kernels

        write_final_code(algorithm_template, final_algorithm, code_file, kernel_file, algorithm.language)
        if algorithm.language == 'F90':
            write_module(self.module, module_file)
            module_file.close()
        code_file.close()
        kernel_file.close()

        latex.write_footer()
        latex.close()

        return

    def sanity_check(self, algorithm):
        """ Perform sanity checks on the system's properties. """

        if len(algorithm.bcs) == self.ndim:
            pass
        elif (len(algorithm.bcs) > self.ndim):
            raise ValueError('There are more boundary conditions than the number of dimensions')
        elif (len(algorithm.bcs) < self.ndim):
            raise ValueError('There are less boundary conditions than the number of dimensions')

        if len(algorithm.explicit_filter) == self.ndim:
            pass
        elif (len(algorithm.explicit_filter) > self.ndim):
            raise ValueError('There are more options for filter than the number of dimensions')
        elif (len(algorithm.explicit_filter) < self.ndim):
            raise ValueError('There are less options for filter than the number of dimensions')

        return

    def compile(self, algorithm, simulation_parameters):
        if algorithm.language == 'OPSC':
            # Translate the generated code using the OPSC translator.
            LOG.debug("Translating OPSC code...")
            try:
                from ops_translator.c import ops
            except ImportError as e:
                LOG.error("Could not import the OPS translator. Aborting...")
                LOG.exception(e)
                return

            try:
                ops.main(["%s.cpp" % simulation_parameters["name"]])
            except Exception as e:
                LOG.error("Something went wrong when calling the OPS translator.")
                LOG.exception(e)
        return

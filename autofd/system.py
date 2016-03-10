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
        self.deltas = et.ndimarray(et.get_indexed_object(len(shape)))
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
        if temporal.scheme == "Forward":
            if temporal.order == 1:
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
    def __init__(self, equations, ranges, computation):
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
        self.computation_type = computation
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
        self.classify_indexed_objects()
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
        coord = et.ndimarray(et.get_indexed_object(len(grid.shape)))
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
            computations.append(ComputationalKernel(eqs, ranges, "Formula Evaluation"))
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
                    computations.append(ComputationalKernel(eq, ranges[number], "Derivative Evaluation"))
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
                    computations.append(ComputationalKernel(eqs, local_range, "Temporary formula Evaluation"))
                    for eq in eqs:
                        newder = der.subs(eq.rhs,eq.lhs)
                    rhs = SD.get_derivativeformula(newder,evals[der].formula[0])
                    eq = Eq(evals[der].work,rhs)
                    computations.append(ComputationalKernel(eq, ranges[number], "Derivative Evaluation"))
            else:
                newder = der
                if all(subev == None for subev in subevals[number]):
                    for req in require[number]:
                        newder = newder.subs(req,evals[req].work)
                else:
                    raise NotImplementedError("Sub evaluations in a mixed derivative")
                rhs = SD.get_derivativeformula(newder,evals[der].formula[0])
                eq = Eq(evals[der].work,rhs)
                computations.append(ComputationalKernel(eq, ranges[number], "Nested Derivative evaluation"))
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
        # The final computations of the residual (change in the rhs terms of the equations)
        # The residual equations are also named as work arrays
        residual_eqs = []
        residual_arrays = []
        for eq in updated_eq:
            wk = grid.work_array('%s%d'%(wkarray,wkind)); wkind = wkind + 1
            residual_arrays.append({eq.lhs:wk})
            residual_eqs.append(Eq(wk,eq.rhs))
        eval_range = [tuple([0,s]) for s in grid.shape]
        computations.append(ComputationalKernel(residual_eqs, eval_range, "Residual of equation"))
        # The residual arrays are tracked to do the computations temporally
        self.computations = computations
        self.residual_arrays = residual_arrays
        return
class TemporalSolution():
    def __init__(self,temporal, grid, const_dt, Spatialsolution):
        if const_dt:
            dt = EinsteinTerm('deltat')
            dt.is_constant = True; dt.is_commutative = True
        else:
            raise NotImplementedError("Varying delta t is not implemented in the code")
        self.nstages = temporal.order
        if temporal.scheme == "Forward" and self.nstages == 1:
            self.coeff = None
            pass
        elif temporal.scheme == "RungeKutta" and self.nstages == 3:
            self.coeff = RungeKutta(self.nstages)
            pass
        else:
            raise ValueError("Only 1st order Forward or RungeKutta 3 order temporal schemes are allowed")
        eqs = []
        # Any computations at the start of the time step generally save equations
        self.start = []
        self.computations = []
        out = []
        for soln in Spatialsolution.residual_arrays:
            out.append(self.time_derivative(soln.keys()[0].args[0], dt,soln[soln.keys()[0]], grid))
        if self.nstages !=1:
            start = [o[-1] for o in out]
            range_ofevaluation = [tuple([0+grid.halos[i][0],s+grid.halos[i][1]]) for i,s in enumerate(grid.shape)]
            self.start.append(ComputationalKernel(start,range_ofevaluation, "Save equations"))
            range_ofevaluation = [tuple([0,s]) for i,s in enumerate(grid.shape)]
            # these are the update equations of the variables at t + k where k is rk loop
            eqs = [o[0] for o in out]
            self.computations.append(ComputationalKernel(eqs,range_ofevaluation, "Rk new (subloop) update"))
            eqs = [o[1] for o in out]
            self.computations.append(ComputationalKernel(eqs,range_ofevaluation, "RK old update"))
        else:
            self.start = None
            range_ofevaluation = [tuple([0,s]) for i,s in enumerate(grid.shape)]
            self.computations.append(ComputationalKernel(out,range_ofevaluation, "RK old update"))
        return
    def time_derivative(self, fn, dt, residual, grid):
        out = fn
        if self.nstages == 1:
            eqn = Eq(fn, fn + dt*residual, evaluate=False)
        elif self.nstages == 3:
            old = grid.work_array('%s_old'%fn.base)
            eqn_fn = Eq(fn, old + self.coeff.new*residual, evaluate=False)
            eqn_old = Eq(old, old + self.coeff.old*residual, evaluate=False)
            saveeq = Eq(old, fn)
            eqn = [eqn_fn,eqn_old, saveeq]
        return eqn
class RungeKutta():
    def __init__(self, order):
        self.stage = Symbol('stage', integer=True)
        self.old = IndexedBase('rkold')
        self.new = IndexedBase('rknew')
        self.old = self.old[self.stage]
        self.new = self.new[self.stage]
        if order == 3:
            self.coeffs = {}
            self.coeffs[self.old] = [-1, 2, -1]
            self.coeffs[self.new] = [-1, 4, -6, 4, -1]
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

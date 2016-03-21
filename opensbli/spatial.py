#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

import os
import subprocess
from sympy import *
from sympy.parsing.sympy_parser import *

import logging
LOG = logging.getLogger(__name__)

from .equations import *
from .opsc import *
from .kernel import *

class SpatialDerivative(object):

    """ The spatial derivatives of an arbitrary function 'F'
    on the numerical grid with the provided spatial scheme.
    
    For a wall boundary condition this will have a dependency on the grid range. """
    
    def __init__(self, spatial_scheme, grid, max_order):
        """ Initialise the spatial derivative, which gives the equations
        of spatial Derivatives for the various combinations of the spatial scheme and order of accuracy.
        
        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :arg int max_order: The maximum order of the derivative in the function.
        :returns: None
        """
        
        # FIXME: The stencil should be formula dependant
        self.stencil = self.create_stencil(spatial_scheme, grid)
        
        self.derivatives = []
        self.derivative_direction = grid.indices
        self.deltas = grid.deltas
        base = IndexedBase('f',shape = grid.shape)
        base.is_grid = True; base.is_constant = False
        fn = base[grid.indices]
        self.fn = fn
        
        self.create_derivative_formulas(fn, max_order, grid)
        return

    def create_stencil(self, spatial_scheme, grid):
        """ Create the computational stencil referenced by the grid indices.
        
        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :returns: A list of solution points indexed on the grid describing the computational stencil.
        :rtype: list
        """
        
        stencil = [[] for dim in grid.shape]
        
        for dim, val in enumerate(grid.shape):
            if spatial_scheme.name == "central":
                points = list(i for i in range(-spatial_scheme.order/2, spatial_scheme.order/2+1)) # The local indices of each point in the stencil (in dimension 'dim').
                grid.halos.append((-spatial_scheme.order/2, spatial_scheme.order/2)) # The range of the indices of the stencil at the boundary which includes the halo points.
            else:
                raise NotImplementedError("Only central difference schemes are supported.")
            stencil[dim] = [grid.indices[dim] + i for i in points] # The indices in the grid offset by the local indices of each stencil point (e.g. i-1, i, i+1)
        return stencil
        
    def create_derivative_formulas(self, fn, max_order, grid):
        """ Create the formulas for the derivatives of the given function,
        based on the stencil pattern provided.
        
        :arg fn: The function whose differential needs computing.
        :arg int max_order: The maximum order of the derivative in the function.
        :arg grid: The numerical grid of solution points.
        :returns: None
        """
    
        derivatives = []
        derivatives += [fn] # FIXME: Later change this to interpolation
        derivative_formula = []
        derivative_formula += [fn]
        kernels = []
        kernels += [fn]
        
        # Create derivatives for each order, up to the maximum order of derivative that is present in the equation.
        for order in range(1, max_order+1):
            shape = tuple([len(grid.indices) for ind in range(order)])
            array = MutableDenseNDimArray.zeros(*shape)
            fdarray = MutableDenseNDimArray.zeros(*shape)
            derivative_kernel = MutableDenseNDimArray.zeros(*shape)
            
            for ind in np.ndindex(*array.shape):
                # Arguments to the derivative in terms of grid indices
                derivative_args = [grid.indices[i] for i in ind]
                # The derivative kernel's name
                name = [str(arg) for arg in ind]
                name = "[%d][%s]" % (order, ','.join(name))
                derivative_kernel[ind] = Symbol(name)
                
                # Find the finite difference formula
                array[ind] = fn.diff(*derivative_args)
                if order == 1 or len(set(derivative_args)) == 1:
                    fdarray[ind] = as_finite_diff(array[ind], self.stencil[ind[0]])*pow(grid.deltas[ind[0]],-order)
                else:
                    new_derivative = array[ind].subs(derivatives[order-1][ind[:-1]], derivative_formula[order-1][ind[:-1]])
                    fdarray[ind] = as_finite_diff(new_derivative, self.stencil[ind[-1]], wrt=grid.indices[ind[-1]])*pow(grid.deltas[ind[-1]],-1)
                    
            derivatives.append(array)
            derivative_formula.append(fdarray)
            kernels.append(derivative_kernel)
            
        self.derivatives = derivatives
        self.derivative_formula = derivative_formula
        self.derivative_kernel = kernels
        return
        
    def get_derivative_formula(self, derivative):
        """ Return the formula for the derivative. For getting a symbolic derivative for a general function, use get_derivative.
        Used for ceval stuff.
        
        :arg derivative: The derivative you want to get the formula for.
        :returns: The derivative's formula.
        """
        
        order = len(derivative.args[1:])
        indices = []
        for arg in derivative.args[1:]:
            indices = indices + [self.derivative_direction.index(arg)]
        if order == 1 or len(set(indices)) == 1:
            formula = as_finite_diff(derivative, self.stencil[indices[0]])*pow(self.deltas[indices[0]], -order)
        else:
            lower_derivative = Derivative(derivative.args[0], *indices[:-1])
            raise ValueError("First update the derivative of %s before calling %s" % (lower_derivative, derivative))
        return formula
        
    def get_derivative(self, derivative):
        """ Return a tuple to which the derivative formula exists in
        the already-evaluated derivatives.
        
        :arg derivative: The derivative you want to get the formula for.
        """
        
        order = len(derivative.args[1:])
        indices = []
        for arg in derivative.args[1:]:
            indices = indices + [self.derivative_direction.index(arg)]
        general_formula = []
        subevals = []
        requires = []
        if order == 1 or len(set(indices)) == 1:
            general_formula += [order,tuple(indices)]
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += [derivative.args[0]]
        else:
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))

            else:
                subevals += [None]
                requires += [derivative.args[0]]
            general_formula += [order-1, tuple([indices[-1]])]
            requires += [self.derivatives[order-1][indices[:-1]].subs(self.fn, derivative.args[0])]

        return general_formula, subevals, requires


class Evaluations(object):

    def __init__(self, lhs, rhs, requires, subevals = None, wk=None):
    
        self.store = True
        
        if isinstance(lhs, Derivative):
            self.is_derivative = True
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
            self.is_derivative = False
            
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


class SpatialDiscretisation(object):

    """ The spatial discretisation using the provided scheme on the provided grid. """

    def __init__(self, equations, formulas, grid, spatial_scheme):
        """ Perform the spatial discretisation. """
    
        all_equations = flatten(list(e.expanded for e in equations))
        all_formulas = flatten(list(e.expanded for e in formulas))
        max_order = maximum_derivative_order(all_equations)
        
        spatial_derivative = SpatialDerivative(spatial_scheme, grid, max_order)
        
        grid_arrays = {}
        range_used = {}
        grid_variables, variable_count = get_indexed_grid_variables(all_equations + all_formulas)
        for atom in grid_variables:
            grid_arrays[atom] = indexed_by_grid(atom, grid)
        
        # Get the SymPy Derivative objects in all of the equations.
        spatial_derivatives, derivative_count, time_derivatives = self.get_spatial_derivatives(all_equations + all_formulas)
        
        # Define the formulas on the grid. This is substituting the old with new.
        # TODO: a sanity check of the formulas, i.e. remove all the formulas that
        # are not used in the equations.
        evals = {}
        for formula in all_formulas:
            out = formula
            for atom in formula.atoms(Indexed):
                out = out.subs(atom, grid_arrays[atom])
            evaluated = Evaluations(out.lhs, out.rhs, list(out.rhs.atoms(Indexed)), None, out.lhs)
            evals[out.lhs] = evaluated
            
        # TODO: a way of passing the coordinates
        term = EinsteinTerm('x_i')
        term.is_constant = True
        coord = term.get_array(term.get_indexed(len(grid.shape)))
        coord = coord.tolist()
        
        # Work array is always named as wk
        work_array = 'wk'; work_array_index = 0;
        for derivative in spatial_derivatives:
            out = derivative # Modify the derivative to be a derivative on grid
            wk = grid.work_array('%s%d'%(work_array,work_array_index)); work_array_index = work_array_index +1
            for atom in derivative.atoms(Indexed):
                out = out.subs(atom, grid_arrays[atom])
            for arg in out.args[1:]:
                out = out.subs(arg, grid.indices[coord.index(arg)])
            general_formula, subevals, requires = spatial_derivative.get_derivative(out)
            grid_arrays[derivative] = out
            evaluated = Evaluations(out,general_formula,requires, subevals, wk)
            evals[out] = evaluated
            
        # We will assume that all the functions in time derivative are known at the start
        known = [grid_arrays[d.args[0]] for d in time_derivatives]
        for val in known:
            evaluated = Evaluations(val,val, None, None,val)
            evals[val] =  evaluated
            
        # Sort the Formulas
        order_of_evaluations = [grid_arrays[d.args[0]] for d in time_derivatives]
        order_of_evaluations = sort_evaluations(order_of_evaluations,evals, Indexed)
        
        # Sort the derivatives
        order_of_evaluations = sort_evaluations(order_of_evaluations,evals, Derivative)
        
        # Update the range of evaluations for each evaluation
        range_of_evaluation(order_of_evaluations, evals,grid, spatial_derivative)
        
        # Now define a Kernel for each of the evaluations

        # All the variables (IndexedObjects) in the equations excluding those which have a time derivative are stored into a kernel
        forms = [ev for ev in order_of_evaluations if isinstance(ev, Indexed) and ev not in known]
        ranges = [evals[ev].evaluation_range for ev in forms]
        subevals = flatten([evals[ev].subevals for ev in forms])
        subeval_truth = [ev == None for ev in subevals]
        
        # Check if all the ranges of evaluation are the same for the formulas.
        range_truth = [ranges[0][i] == val[i] for val in ranges for i in range(len(ranges[0]))]
        computations = []
        eqs = []
        eqs = [Eq(evals[ev].work, evals[ev].formula) for ev in forms]
        if forms:
            # If same range then combine them into a single computation else store into different computations
            if all(range_truth) and all(subeval_truth):
                computations.append(Kernel(eqs, ranges[0], "Formula Evaluation"))
            else:
                for number,eq in enumerate(eqs):
                    computations.append(Kernel(eq, ranges[number]))
                    
        # Now process the Derivatives
        # TODO: This can be moved out into a seperate function. Which can be used for Diagnostics/Generalised coordinate equations evaluations
        derivatives = [ev for ev in order_of_evaluations if isinstance(ev, Derivative) and ev not in known]
        ranges = [evals[ev].evaluation_range for ev in derivatives]
        subevals = [evals[ev].subevals for ev in derivatives]
        require = [evals[ev].requires for ev in derivatives]
        for number, derivative in enumerate(derivatives):
            if not any(isinstance(req, Derivative) for req in require[number]):
                if all(subev == None for subev in subevals[number]):
                    rhs = spatial_derivative.get_derivative_formula(derivative)
                    eq = Eq(evals[derivative].work,rhs)
                    computations.append(Kernel(eq, ranges[number], "Derivative Evaluation"))
                else:
                    # Store into temporary array the sub evaluation
                    eqs = []
                    temp_work_array_index = work_array_index
                    for subev in subevals[number]:
                        wk = grid.work_array('%s%d' % (work_array, temp_work_array_index))
                        temp_work_array_index += 1
                        for req in require[number]:
                            local_range = evals[req].evaluation_range
                            subev = subev.subs(req, evals[req].work)
                        eqs.append(Eq(wk,subev))
                    computations.append(Kernel(eqs, local_range, "Temporary formula Evaluation"))
                    for eq in eqs:
                        new_derivative = derivative.subs(eq.rhs, eq.lhs)
                    rhs = spatial_derivative.get_derivative_formula(new_derivative)
                    eq = Eq(evals[derivative].work,rhs)
                    computations.append(Kernel(eq, ranges[number], "Derivative Evaluation"))
            else:
                new_derivative = derivative
                if all(subev == None for subev in subevals[number]):
                    for req in require[number]:
                        new_derivative = new_derivative.subs(req,evals[req].work)
                else:
                    raise NotImplementedError("Sub evaluations in a mixed derivative")
                rhs = spatial_derivative.get_derivative_formula(new_derivative)
                eq = Eq(evals[derivative].work,rhs)
                computations.append(Kernel(eq, ranges[number], "Nested Derivative evaluation"))

        # All the spatial computations are evaluated by this point now get the updated equations
        updated_equations = [e for e in all_equations]
        for equation_number, equation in enumerate(updated_equations):
            spatial_derivatives, derivative_count, time_derivatives = self.get_spatial_derivatives([equation])
            grid_variables, variable_count = get_indexed_grid_variables([equation])
            spatial_derivatives = (sorted(spatial_derivatives, cmp = decreasing_order))
            # substitute spatial derivatives first
            for var in spatial_derivatives + grid_variables:
                new = evals[grid_arrays[var]].work
                updated_equations[equation_number] = updated_equations[equation_number].subs(var, new)
                
        # The final computations of the residual (change in the RHS terms of the equations).
        # The residual equations are also named as work arrays.
        # The residual arrays are tracked for use in the evaluation of the temporal scheme.
        residual_equations = []
        residual_arrays = []
        for e in updated_equations:
            work_array = grid.work_array('%s%d' % (work_array, work_array_index))
            work_array_index += 1
            residual_arrays.append({e.lhs : work_array})
            residual_equations.append(Eq(work_array, e.rhs))
        evaluation_range = [tuple([0, s]) for s in grid.shape]
        computations.append(Kernel(residual_equations, evaluation_range, "Residual of equation"))

        # Update the required computations and residual arrays
        self.computations = computations
        self.residual_arrays = residual_arrays
        return

    def get_spatial_derivatives(self, equations):
        """ Return all the spatial Derivative terms in the equations.
        Any equations involving Derivative objects in terms of the time 't' are handled separately.
        
        :arg equations: A list of equations to search.
        :returns: All of the spatial Derivative objects, the count of the number of times a particular Derivative object is used, and all of the temporal Derivative objects.
        """
    
        derivatives = []
        count = {}
        time_derivatives = []
        
        for eq in equations:
            pot = preorder_traversal(eq)

            for p in pot:
                if p in derivatives:
                    pot.skip()
                    count[p] = count[p]+1
                    continue
                elif isinstance(p, Derivative):
                    if all(arg != EinsteinTerm('t') for arg in p.args):
                        pot.skip()
                        derivatives.append(p)
                        count[p] = 1
                    else:
                        pot.skip()
                        time_derivatives.append(p)
                else:
                    continue
                    
        return derivatives, count, time_derivatives
    


class GridBasedInitialisation(object):

    """ Initialise the equations on the grid of solution points.  """

    def __init__(self, grid, ics):
        """ Create the initialisation kernels.
        
        :arg grid: The numerical grid of solution points.
        :arg list ics: A list of initial condition formulas.
        :returns: None
        """
    
        self.computations = []
        initialisation_equation =[]
        for ic in ics:
            initialisation_equation.append(parse_expr(ic, local_dict = {'grid':grid}))
        range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
        
        self.computations.append(Kernel(initialisation_equation, range_of_evaluation, "Initialisation"))
        
        return
        

def range_of_evaluation(order_of_evaluations, evaluations, grid, sdclass):
    """ First the ranges of derivatives are updated, then other ranges are updated. """
    
    derivatives = []
    for ev in order_of_evaluations:
        if isinstance(ev, Derivative):
            derivatives.append(ev)
        evaluations[ev].evaluation_range = [tuple([0, s]) for s in grid.shape]
        
    # Update the range for the derivatives
    grouped_derivatives = group_derivatives(derivatives)
    for key, value in grouped_derivatives.iteritems():
        for val in value:
            require = evaluations[val].requires
            formula  = evaluations[val].formula
            direction = formula[1][0]
            halos = grid.halos[direction]
            for req in require:
                erange = list(evaluations[req].evaluation_range[direction])
                if erange[0] == 0 and erange[1] == grid.shape[direction]:
                    erange[0] = erange[0]+halos[0]
                    erange[1] = erange[1]+halos[1]
                evaluations[req].evaluation_range[direction] = tuple(erange)
                
    # Update the range for the formulas
    for ev in order_of_evaluations:
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


def group_derivatives(derivatives):
    """ For each SymPy Derivative term, get the term that it is operating on and
    formula (function, Derivative(function)) pairs in a dictionary.
    
    :arg list derivatives: A list of all the Derivative terms to consider.
    :returns: A dictionary of the (function, Derivative(function)) pairs.
    :rtype: dict
    """

    derivative_dict = {}
    
    for derivative in derivatives:
        if derivative.args[0] in derivative_dict.keys():
            derivative_dict[derivative.args[0]] += [derivative]
        else:
            derivative_dict[derivative.args[0]] = [derivative]

    for key,value in derivative_dict.iteritems():
        if len(value) > 1:
            derivative_dict[key] = (sorted(value, cmp=increasing_order))

    return derivative_dict


def indexed_by_grid(variable, grid):
    """ Convert a variable/function or Indexed object to an Indexed object indexed by the Grid indices.
    
    :arg variable: The variable to convert to a Grid-based Indexed variable
    :arg grid: The numerical Grid of solution points.
    :returns: An Indexed variable, which is the same variable as the one provided, but is indexed by the Grid indices.
    :rtype: sympy.Indexed
    """
    
    if isinstance(variable, Indexed):
        base = IndexedBase('%s' % variable.base)
    elif isinstance(variable, Function):
        base = IndexedBase('%s' % variable.func)
    else:
        raise ValueError("Only functions or Indexed Objects are supported", variable)
    base.is_grid = True; base.is_constant = False
    return base[grid.indices]


def get_indexed_grid_variables(equations):
    """ Return all the variables in the equations that are Indexed on the Grid.
    
    :arg list equations: A list of SymPy equations to consider.
    :returns: A list of the variables in the equations that are Indexed on the Grid, and also the count of all the specific terms in the equations (Indexed or not).
    :rtype: (list, int)
    """

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

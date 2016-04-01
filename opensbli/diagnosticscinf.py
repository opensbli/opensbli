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

'''
For writing the variables that are already evaluated we are using FileIO
The various diagnostics to be supported are
1. Compute some variables and write to files (Eg. Statistics, vorticity, etc etc)
    this is sub divided into
    a. Do any computations if exist
    b. These can be after n time steps or at the end of the simulation
        i. Check point b, if after n iterations write if statement
        ii. Else at the end of simulation write it

2. Reduction of a variable i.e already evaluated or evaluation and write to a file
    this can be sub divded into
    a. Computation
    b. Reduction operation
    c. Writing to file (create a new file or append to a file)

Now the inputs to the class,
These are the equations expanded using Einstein expansion
type of Diagnostics time diagnostics, time average, reduction after n time steps,

'''

from sympy.tensor import IndexedBase, Indexed
from .spatial import *
from .opsc import *
class Diagnostics(object):
    computations = None
    after_time_computations = None
    reduction_type = None
    save_after = 0
    file_name = ''
class ReductionVariable(Symbol):
    def __new__(self,var):
        self = Symbol.__xnew__(self, var)
        return self
class Reduction():

    def __init__(self, grid, equations, formulas, spatial_scheme):
        '''
        First try out the reduction code for finding Kinetic energy
        eq will be 0.5*(u_j**2). as of now only the equations
        evaluated_variables should be the ones that are known at the starting of diagnostics
        these will be the prognostic variables
        fname will be automatic??


        '''
        #self.type = reduction_type
        self.computations = []
        all_equations = flatten(equations)
        all_formulas = flatten(formulas)
        #self.fname = fname
        #self.save_after = after
        all_formulas = get_used_formulas(all_formulas, all_equations)
        spatial_derivatives, time_derivatives = get_derivatives(all_equations + all_formulas)
        evaluations = {}
        spatial_derivative = SymDerivative(spatial_scheme, grid)
        evaluations = create_formula_evaluations(all_formulas, evaluations)
        evaluations = create_derivative_evaluations(spatial_derivatives,evaluations, spatial_derivative)
        pprint(evaluations)
        pprint(all_equations)
        ke = ReductionVariable("ke")
        pprint(type(ke))
        return
    def sum_reduction(self, equations):
        """

        """

        return

class SymDerivative(object):
    """ The spatial derivatives of an arbitrary function 'F'
    on the numerical grid with the provided spatial scheme.

    For a wall boundary condition this will have a dependency on the grid range. """

    def __init__(self, spatial_scheme, grid):
        """ Initialise the spatial derivative, which gives the equations
        of spatial Derivatives for the various combinations of the spatial scheme and order of accuracy.

        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :returns: None
        """
        self.derivative_direction = grid.mapedindices
        self.deltas = grid.deltas
        self.points = spatial_scheme.points
        return
    def get_derivative_formula(self, derivative):
        """
        This returns the formula for the derivative function
        """
        order = len(derivative.args[1:])
        indices = [self.derivative_direction.index(arg) for arg in derivative.args[1:]]
        if order == 1 or len(set(indices)) == 1:
            wrt = indices[0]
            gridpoints = [wrt + i for i in self.points]
            formula = apply_finite_diff(order, gridpoints, [derivative.expr.subs({wrt: x}) for x in gridpoints], wrt)
        elif order == 2:
            # Do the derivative of derivative
            raise NotImplementedError("Derivatives of order == 2  Mixed is not implemented")
        else:
            raise NotImplementedError("Derivatives of order > 2 are not implemented")
        return
    def get_derivative(self, derivative):
        """ Return a tuple to which the derivative formula exists in
        the already-evaluated derivatives.

        :arg derivative: The derivative you want to get the formula for.
        """
        pprint(derivative)
        order = len(derivative.args[1:])
        indices = []
        for arg in derivative.args[1:]:
            pprint(arg)
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
                requires += [derivative.args[0].atoms(Indexed)]
        else:
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += [derivative.args[0]]
            general_formula += [order-1, tuple([indices[-1]])]
            requires += [Derivative(derivative.args[0],derivative.args[1:-1])]

        return general_formula, subevals, requires

'''
Here I rewrite the stuff in spatial.py, esentially dividing the big spatial.py into smaller stuff
If working fine these should be used back in spatial.py

All the definitions here will be moved to a seperate file that can be used both in
spatial and other places

the essense would be to
1. Get all the formulas used in the equations
2. Get the derivatives used in the equations
'''
def convert_equations_to_grid(equations, grid):
    """
    Converts the indexed objects in the equations to indexed objects on the grid
    # NOW THIS IS NOT REQUIRED
    """
    grid_equations = []

    for eq in equations:
        substitutions = {}
        variables, count = get_indexed_grid_variables([eq])
        for var in variables:
            substitutions[var] = indexed_by_grid(var, grid)
        grid_equations.append(eq.subs(substitutions))
    pprint(grid_equations)
    return grid_equations

def get_used_formulas(formulas, equations):
    '''
    This returns the formulas used in the equations.
    '''
    variables, count = get_indexed_grid_variables(equations)

    formulas = dict(zip([form.lhs for form in formulas], [form.rhs for form in formulas]))

    used_formulas = [Eq(var,formulas[var]) for var in variables if var in formulas.keys()]

    return used_formulas
def create_derivative_evaluations(spatial_derivatives, evals, symderivative):
    """
    Derivative computations are evaluated seperately as they sometimes require evaluation of
    temporary work arrays
    """
    for out in spatial_derivatives:
        general_formula, subevals, requires = symderivative.get_derivative(out)
        evaluated = Evaluations(out, general_formula, requires, subevals, None)
        evals[out] = evaluated
    return evals

def create_formula_evaluations(all_formulas, evals):
    """
    Creates computation class for all the equations provided.
    This is used for creating any computation (Diagnostics, Spatial, Temporal, Shock-capturing)
    """
    for formula in all_formulas:
        evaluated = Evaluations(formula.lhs, formula.rhs, list(formula.rhs.atoms(Indexed)), None, None)
        evals[formula.lhs] = evaluated

    return evals

'''
Simlar to the classes create stuff
'''

#def indexed_by_grid(variable, grid):
    #""" Convert a variable/function or Indexed object to an Indexed object indexed by the Grid indices.

    #:arg variable: The variable to convert to a Grid-based Indexed variable
    #:arg grid: The numerical Grid of solution points.
    #:returns: An Indexed variable, which is the same variable as the one provided, but is indexed by the Grid indices.
    #:rtype: sympy.Indexed
    #"""

    #if isinstance(variable, Indexed):
        #base = IndexedBase('%s' % variable.base)
    #elif isinstance(variable, Function):
        #base = IndexedBase('%s' % variable.func)
    #else:
        #raise ValueError("Only functions or Indexed Objects are supported", variable)
    #base.is_grid = True; base.is_constant = False
    #return base[grid.indices]

def sort_evaluations(order, evaluations, typef):
    """ Sort the evaluations based on the requirements of each term. For example, if we have
    the primitive variables p, u0, u1, and T, then the pressure p may depend on the velocity u0 and u1, and T may depend on p,
    so we need this be evaluate in the following order: u0, u1, p, T.

    :arg list order: The list of terms to sort.
    :arg evaluations: The evaluation information, containing dependency information.
    :arg typef: The type of term to sort.
    :returns: A list of ordered terms.
    :rtype: list
    """

    for key in evaluations.keys():
        if isinstance(key, typef) and not key in order:
            if all(ev in order for ev in evaluations[key].requires):
                order.append(key)
            else:
                for val in evaluations[key].requires:
                    if not val in order:
                        sort_evaluations(order, {val:evaluations[val]}, typef)
                order.append(key)
    return order
#def set_range_of_evaluations(order_of_evaluations, evaluations, grid):
    #""" Set the evaluation ranges of each Evaluation object based on the shape of the grid of points.
    #First the ranges of derivatives are updated, then other ranges are updated. """

    #derivatives = []
    #for ev in order_of_evaluations:
        #if isinstance(ev, Derivative):
            #derivatives.append(ev)
        #evaluations[ev].evaluation_range = [tuple([0, s]) for s in grid.shape]

    ## Update the range for the derivatives
    #grouped_derivatives = group_derivatives(derivatives)
    #for key, value in grouped_derivatives.iteritems():
        #for val in value:
            #require = evaluations[val].requires
            #formula  = evaluations[val].formula
            #direction = formula[1][0]
            #halos = grid.halos[direction]
            #for req in require:
                #erange = list(evaluations[req].evaluation_range[direction])
                #if erange[0] == 0 and erange[1] == grid.shape[direction]:
                    #erange[0] += halos[0]
                    #erange[1] += halos[1]
                #evaluations[req].evaluation_range[direction] = tuple(erange)

    ## Update the range for the formulas this might require some modifications
    #for ev in order_of_evaluations:
        #if isinstance(ev, Indexed):
            #require = evaluations[ev].requires
            #if require:
                #for req in require:
                    #evaluations[req].evaluation_range = evaluations[ev].evaluation_range

    #return

def get_derivatives(equations):
    """ Return all the spatial Derivative terms in the equations.
    Any equations involving Derivative objects in terms of the time 't' are handled separately.

    :arg equations: A list of equations to search.
    :returns: All of the spatial Derivative objects and all of the temporal Derivative objects.
    """

    derivatives = []
    time_derivatives = []

    for eq in equations:
        pot = preorder_traversal(eq)

        for p in pot:
            if p in derivatives:
                pot.skip()
                continue
            elif isinstance(p, Derivative):
                if all(arg != EinsteinTerm('t') for arg in p.args):
                    pot.skip()
                    derivatives.append(p)
                else:
                    pot.skip()
                    time_derivatives.append(p)
            else:
                continue

    return derivatives, time_derivatives

#def get_indexed_grid_variables(equations):
    #""" Return all the variables in the equations that are Indexed on the Grid.

    #:arg list equations: A list of SymPy equations to consider.
    #:returns: A list of the variables in the equations that are Indexed on the Grid, and also the count of all the specific terms in the equations (Indexed or not).
    #:rtype: (list, int)
    #"""

    #variables = []
    #count = {}
    #for eq in equations:
        #pot = preorder_traversal(eq)
        #for p in pot:
            #if p in variables:
                #pot.skip()
                #count[p] = count[p]+1
                #continue
            #elif isinstance(p, Indexed):
                #pot.skip()
                #variables.append(p)
                #count[p] = 1
            #else:
                #continue
    #return variables, count
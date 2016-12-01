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
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>

from sympy import *

from .opensbliobjects import CoordinateObject, CoordinateObject, DataObject, DataSet
from .opensblifunctions import *
classdict = {Symbol('Central'): CentralDerivative}

# from .latex import LatexWriter
from sympy.tensor.array import MutableDenseNDimArray
import sys


class Transformation(object):
    def __init__(self, ndim, coordinate, parameters):
        self.ndim = ndim
        self.stretching_metric = [param[0] for param in parameters] # Get whether true or false fr streching
        self.curvilinear_metric = [param[1] for param in parameters] # Get T/F for curvilinear
        cart = CoordinateObject('%s_i'%(coordinate))
        curv = CoordinateObject('%s_i'%('xi'))
        cartesian_coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        #cartesian_coordinates = [CoordinateObject('%s%d'%(coordinate, ind)) for ind in range(ndim)]
        curvilinear_coordinates = [curv.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        self.cart_directions = range(self.ndim)
        self.curvilinear_coordinates = curvilinear_coordinates
        self.cartesian_coordinates = cartesian_coordinates
        self.cart_to_curvilinear_functions = Matrix(ndim, ndim,lambda i,j:curvilinear_coordinates)
        self.curvilinear_to_cart_functions = Matrix(ndim, ndim,lambda i,j:cartesian_coordinates)
        self.update_curvilinear_function(self.curvilinear_metric)
        self.update_stretching_function(self.stretching_metric)
        self.update_functions()
        """
        Is there an easy way to parse them as strings and substitute zeros??
        First Derivatives
        jacobians = CD(xi_i,x_j)
        FDs = CD(f(xi_i(x0,x1
        """

        self.evaluate_transform_matrix()
        # Generate transform indexing
        #self.generate_transformation_indices()
        self.evaluate_metrics_fd()
        # Generate derivatives with respect to each of the spatial dimensions
        A1, B1 = self.generate_first_derivatives()
        self.generate_second_derivatives(A1,B1)
        # First derivative object
        DataSet.dimensions = ndim
        self.coordinateobjectsmapping()
        self.evaluations = {}
        order = 1
        scheme =  Symbol('Central')
        convertscheme = lambda x:self.convert_to_scheme(scheme,x)
        self.evaluations[order] = TransformationEvaluation(order)
        pprint(self.FD)
        pprint(self.adjointJ)
        pprint(self.detJ)
        exit()
        self.evaluations[order].transformed_derivative = self.FD.applyfunc(convertscheme)
        # Convert the scheme to the type of discretisation scheme
        #pprint(self.adjointJ.applyfunc(convertscheme))
        #for t in self.evaluations[order].transformed_derivative:
            #for d in (t.atoms(CentralDerivative)):
                #pprint(d.__dict__)
        #pprint(self.detJ)
        self.evaluations[order].symbols = self.FD_symbols.applyfunc(convertscheme)
        print "here"
        pprint(self.FD_eval_symbols)
        self.evaluations[order].eval_symbols = self.FD_eval_symbols.applyfunc(convertscheme)
        self.evaluations[order].metric_evaluations = self.adjointJ.applyfunc(convertscheme)
        self.evaluations[order].detJ = self.convert_to_scheme(scheme, self.detJ)
        order =2
        # Second derivative object
        self.evaluations[order] = TransformationEvaluation(order)
        self.evaluations[order].transformed_derivative = self.SD.applyfunc(convertscheme)
        self.evaluations[order].symbols = self.SD_symbols.applyfunc(convertscheme)
        self.evaluations[order].eval_symbols = self.SD_eval_symbols.applyfunc(convertscheme)
        self.evaluations[order].metric_evaluations = self.adjointJ.applyfunc(convertscheme)
        self.evaluations[order].detJ = self.convert_to_scheme(scheme, self.detJ)
        self.create_transform_equations([1,2])

        # self.print_latex(self.SD)
        # pprint(self.SD)
        # # pprint(self.adjointJ)
        return

    def print_latex(self, item):
        for eq in item:
            print latex(eq)

        return

    def update_stretching_function(self, stretch):
        for ind, val in enumerate(stretch):
            if val:
                list1 = self.cart_to_curvilinear_functions[ind,ind]
                list1[ind] = self.curvilinear_coordinates[ind]
                list2 = self.curvilinear_to_cart_functions[ind,ind]
                list2[ind] = self.cartesian_coordinates[ind]
            else:
                list1 = self.cart_to_curvilinear_functions[ind,ind]
                list1[ind] = False
                list2 = self.curvilinear_to_cart_functions[ind,ind]
                list2[ind] = False
        return
    def update_curvilinear_function(self, curv):
        for ind, val in enumerate(curv):
            if not val:
                for d in range(self.ndim):
                    list1 = self.cart_to_curvilinear_functions[ind,d]
                    list2 = self.cart_to_curvilinear_functions[d,ind]
                    list3 = self.curvilinear_to_cart_functions[ind,d]
                    list4 = self.curvilinear_to_cart_functions[d,ind]
                    if d==ind:
                        for dim in range(self.ndim):
                            list1[dim] = False
                            list2[dim] = False
                            list3[dim] = False
                            list4[dim] = False
                    else:
                        for dim in range(self.ndim):
                            list1[dim] = False
                            list2[dim] = False
                            list3[dim] = False
                            list4[dim] = False
                    for d1 in range(self.ndim):
                        list1 = self.cart_to_curvilinear_functions[d,d1]
                        list1[ind] = False
                        list3 = self.curvilinear_to_cart_functions[d,d1]
                        list3[ind] = False
        return
    def update_functions(self):
        for d in range(self.ndim):
            for d1 in range(self.ndim):
                list1 = self.cart_to_curvilinear_functions[d,d1]
                if not any(list1):
                    if d==d1:
                        self.cart_to_curvilinear_functions[d,d1] = self.curvilinear_coordinates[d]
                    else:
                        self.cart_to_curvilinear_functions[d,d1] = 0
                else:
                    inds = [ind for ind,val in enumerate(list1) if not val]
                    list1 = [v for i, v in enumerate(list1) if i not in inds]
                    self.cart_to_curvilinear_functions[d,d1] = self.cartesian_coordinates[d](*list1)
                #  curvilinear to cartesian
                list1 = self.curvilinear_to_cart_functions[d,d1]
                if not any(list1):
                    if d==d1:
                        self.curvilinear_to_cart_functions[d,d1] = self.cartesian_coordinates[d]
                    else:
                        self.curvilinear_to_cart_functions[d,d1] = 0
                else:
                    inds = [ind for ind,val in enumerate(list1) if not val]
                    list1 = [v for i, v in enumerate(list1) if i not in inds]
                    self.curvilinear_to_cart_functions[d,d1] = self.curvilinear_coordinates[d](*list1)

        return
    def evaluate_transform_matrix(self):
        self.cart_to_curv_mat = Matrix(self.ndim,self.ndim,lambda i,j:self.cart_to_curvilinear_functions[i,j].diff(self.curvilinear_coordinates[j]))
        return


    def generate_transformation_indices(self):
        Einstein = lambda x:x.subs({a:CoordinateObject(str(a)) for a in x.atoms(Symbol)})
        remove = lambda x:self.remove_fn(x)
        transformation_indices = self.curvilinear_to_cart_functions.applyfunc(remove)
        transformation_indices = transformation_indices.applyfunc(Einstein)
        shape = transformation_indices.shape
        self.transformation_indices = tuple([transformation_indices[i,j] for i in range(shape[0]) for j in range(shape[1]) if i==j])
        return

    def evaluate_metrics_fd(self):
        """ Function to generate the determinant of the Jacobian and its adjoint to be used by the
        transformations, both are stored to self.
        :args: None
        :returns: None
        """
        self.adjointJ = self.cart_to_curv_mat.adjugate()
        self.detJ = self.cart_to_curv_mat.det()
        # Change function arguments to be consistent with OpenSBLI
        for fn in (self.detJ.atoms(Function)):
            self.detJ = self.detJ.subs(fn,(DataObject(str(fn.func))))
        for i, entry in enumerate(self.adjointJ):
            for fn in entry.atoms(Function):
                self.adjointJ[i] = self.adjointJ[i].subs(fn,(DataObject(str(fn.func))))
        return

    def generate_first_derivatives(self):
        """Create the first derivative and D_xi_x symbol matrices, and the first derivative formula in the new coordinates
        for a generic function f. The values stored in A1, B1 are used to calculate second derivatives in generate_second_derivatives
        , the first derivative components are stored to self.
        :args: None
        :returns: A1: 3x3 Matrix of placeholder symbols D_xi_x, to be passed to the second derivative generation.
        :returns: B1: 3x3 Matrix of corresponding derivatives for the substitutions stored in A1.
        """
        A1 = zeros(self.ndim, self.ndim)
        B1 = zeros(self.ndim, self.ndim)

        first_der_vec = zeros(self.ndim,1)
        # Loop over each Cartesian direction
        for direction in self.cart_directions:
            # Differentiation with respect to x_0, x_1 or x_2
            wrt = self.cartesian_coordinates[direction]
            f = Function('f')
            # Store the wild function as we use it later in transformations
            self.wildfunc = f
            mat_column = self.curvilinear_to_cart_functions[:,direction]
            # Differentiate the functions
            derivatives = f(*mat_column).diff(wrt)
            first_der_vec[direction, 0] = derivatives
            # Populate the substitution matrices
            for j, xi in enumerate(mat_column):
                # Create Jacobian function J
                f_jac = Function('D_xi%dx%d' % (direction, j))
                A1[direction,j] = f_jac(*mat_column)

        self.store_first_ders = first_der_vec
        B1 = Matrix(self.ndim,self.ndim,lambda i,j:self.curvilinear_to_cart_functions[i,j].diff(self.cartesian_coordinates[j]))

        # Remove zero entries in A1 and B1 matrices
        for i in range(A1.shape[0]):
            for j in range(A1.shape[1]):
                if not isinstance(B1[i,j] ,Derivative):
                    A1[i,j] = B1[i,j]
        # Create substitution lists to restore indices and apply to the first derivative formulae
        indexed_FD_subs_list = {}
        FD_subs_list = dict([(b,a) for (b,a) in zip(list(B1), list(A1)) if isinstance(b, Derivative)])

        for key, value in FD_subs_list.iteritems():
            val = (DataObject(str(self.remove_fn(self.remove_fn(value)))))
            indexed_FD_subs_list[key] = val
        indexed_subs_D = lambda x:x.subs(indexed_FD_subs_list)
        first_der_vec = first_der_vec.applyfunc(indexed_subs_D)
        first_der_vec = first_der_vec.applyfunc(self.remove_fn).applyfunc(self.remove_fn)

        # Store the first derivative formulae, D_xi_x symbols and their corresponding derivatives to self
        self.FD = first_der_vec
        self.FD_symbols = A1.applyfunc(self.remove_fn).applyfunc(self.remove_fn)
        self.FD_simple_symbols = A1.applyfunc(self.remove_fn).applyfunc(self.remove_fn)
        self.FD_eval_symbols = B1
        # Convert the B1 matrix saved to FD_symbols to indexed form
        # This is different to the B1 returned by the function to be used by
        # Generate second derivatives
        for i,val in enumerate(list(self.FD_symbols)):
            if isinstance(val,Symbol):
                self.FD_symbols[i] = DataObject(str(val))
        # Return A1, B1 to be used by generate_second_derivatives
        return A1, B1

    def generate_second_derivatives(self, A1, B1):
        """Create the second derivative and SD_xi_x_xi symbol matrices, and the second derivative formulae in the new coordinates
        for a generic function f. All values are stored to self.
        :args: A1: 3x3 Matrix of placeholder symbols D_xi_x, used to calculate the second derivative formulae
        :args: B1: 3x3 Matrix of the corresponding derivatives for the substitutions stored in A1.
        :returns: None
        """
        first_der_vec = self.store_first_ders
        # Create substitution list for the first derivatives and their symbols from A1, B1
        FD_subs_list = dict([(b,a) for (b,a) in zip(list(B1), list(A1)) if isinstance(b, Derivative) ])
        pprint(first_der_vec)
        subs_D = lambda x:x.subs(FD_subs_list)
        second_sder_fn = first_der_vec.applyfunc(subs_D)
        pprint(second_sder_fn)
        second_der_matrix = Matrix(self.ndim,self.ndim, lambda i,j: second_sder_fn[i].diff(self.cartesian_coordinates[j]))
        second_der_matrix = second_der_matrix.applyfunc(self.remove_fn).applyfunc(self.remove_fn)
        pprint(second_der_matrix)
        exit()
        A2 = MutableDenseNDimArray.zeros(self.ndim, self.ndim, self.ndim)
        B2 = MutableDenseNDimArray.zeros(self.ndim, self.ndim, self.ndim)
        # Generate substitution matrices for second derivatives
        for i in range(self.ndim):
            for j in range(self.ndim):
                for k in range(self.ndim):
                    f_res = Symbol('SD_xi%dx%dxi%d' % (i,j,k))
                    A2[i,j,k] = f_res
                    wrt = self.curvilinear_to_cart_functions[k,k]
                    if wrt != 0:
                        B2[i,j,k] = A1[i,j].diff(wrt)

        B2 = B2.applyfunc(self.remove_fn).applyfunc(self.remove_fn)
        non_zero_ders = [(index, b) for (index,b) in enumerate(B2) if isinstance(b, Derivative)]
        zero_indices = [i for (i,b) in non_zero_ders]
        for i,val in enumerate(list(A2)):
            if i not in zero_indices:
                A2[i] = 0
            else:
                A2[i] = DataObject(str(A2[i]))

        new_B = [str(self.remove_fn(self.remove_fn(b))) for a,b in FD_subs_list.items() ]
        subs3 = {CoordinateObject(b):DataObject(str(b)) for b in new_B}
        SD_subs_list = dict(dict([(b,a) for (b,a) in zip(list(B2), list(A2)) if isinstance(b, Derivative)]).items())

        subs2 = {self.remove_fn(a):DataObject(str(self.remove_fn(self.remove_fn(b)))) for a,b in FD_subs_list.items()}

        subs_SD = lambda x:x.subs(SD_subs_list)
        second_der_matrix = second_der_matrix.applyfunc(subs_SD)
        subs_SD = lambda x:x.subs(subs2)

        second_der_matrix = second_der_matrix.applyfunc(subs_SD)
        third_subs = lambda x:x.subs(subs3)
        second_der_matrix = second_der_matrix.applyfunc(third_subs)

        indexed_SD_subs_list = {}
        for key, value in FD_subs_list.iteritems():
            val = (DataObject(str(self.remove_fn(self.remove_fn(value)))))
            indexed_SD_subs_list[key] = val

        indexed_subs_SD = lambda x:x.subs(indexed_SD_subs_list)
        second_der_matrix = second_der_matrix.applyfunc(indexed_subs_SD)

        second_der_matrix = second_der_matrix.applyfunc(self.remove_fn).applyfunc(self.remove_fn)
        self.SD = second_der_matrix
        self.SD_symbols = A2

        for no, val in enumerate((B2)):
            if val.atoms(Derivative):
                val = val.subs({der.args[0]:DataObject(str(der.args[0]))\
                    for der in val.atoms(Derivative) if not isinstance(der, DataObject)})
                B2[no] = val
        self.SD_eval_symbols = B2
        return

    def remove_fn(self,expr):
        #expr = self.apply_Subs(expr)
        for ind in expr.atoms(Function):
            if not isinstance(ind, self.wildfunc):
                expr = expr.subs({ind:CoordinateObject(str(ind.func))})
            else:
                expr = expr.subs({ind:DataObject(str(ind.func))})
        return expr

    def apply_Subs(self,expression):
        substis = list(expression.atoms(Subs))
        for work in substis:
            test = list(zip(work.variables, work.point))
            args = work.expr.args
            newargs = [arg.subs(test) for arg in args]
            work1 = type(work.expr)(*newargs)
            expression = expression.subs(work, work1)
        return expression

    def convert_to_scheme(self, scheme, expr):
        pot = postorder_traversal(expr)
        substitutions = {}
        modifications = []
        for d in pot:
            if isinstance(d, Derivative):
                if d.args[0].atoms(Derivative):
                    modifications += [d]
                else:
                    substitutions[d] = classdict[scheme](d.expr,*d.variables)
            else:
                continue
        expr = expr.xreplace(substitutions)
        modifications = [m.xreplace(substitutions) for m in modifications]
        substitutions = {}
        for m in modifications:
            if not (m.args[0].atoms(Derivative).difference(m.args[0].atoms(classdict[scheme]))):
                expr = expr.xreplace({m:classdict[scheme](m.expr,*m.variables)})
            else:
                raise ValueError("Require better implementation in set schemes")
        # Now convert DataObject to DataSets
        for at in expr.atoms(DataObject):
            obj = DataSet(str(at)) # By default the location of the dataset is node (0,0,0)
            expr = expr.replace(at, obj)
        # set the direction for the coordinate object TODO
        pprint(expr)
        for at in expr.atoms(CoordinateObject):
            #pprint(at)
            at.direction = self.coordinateobjectmap.index(at)

        return expr
    def coordinateobjectsmapping(self):
        remove = lambda x:self.remove_fn(x)
        tranformation_indices = self.curvilinear_to_cart_functions.applyfunc(remove)
        shape = tranformation_indices.shape
        transformation_indices = ([CoordinateObject(str(tranformation_indices[i,j])) for i in range(shape[0]) for j in range(shape[1]) if i==j])
        pprint(transformation_indices)
        self.coordinateobjectmap = transformation_indices
        return

    """ From HERE On these are specific to transformation of equations"""
    def map_transformation_indices(self):
        """ Maps the transformation indices i.e (x0,x1,x2) --> (xi0,xi1,xi2) depending on the input metric
        parameters"""
        remove = lambda x:self.remove_fn(x)
        tranformation_indices = self.curvilinear_to_cart_functions.applyfunc(remove)
        self.tranformation_indices = tranformation_indices
        shape = tranformation_indices.shape
        tranformation_indices = ([CoordinateObject(str(tranformation_indices[i,j])) for i in range(shape[0]) for j in range(shape[1]) if i==j])
        self.tranformation_indices = tranformation_indices
        grid_indices = list(self.indices)
        self.cartesian_to_curvilinear_mapping = dict(zip(self.coordinates , tranformation_indices))
        self.mapped_indices = dict(zip(tranformation_indices+grid_indices , grid_indices + tranformation_indices))
        return

    def transform_first_derivative_vector(self, vector, direction):
        order = 1
        if order in self.evaluations.keys():
            if not self.evaluations[order].used:
                self.evaluations[order].used = True
                self.evaluations[order].Jacobian_indexed = DataObject('Jac', *self.tranformation_indices)
                self.evaluations[order].Jacobian_indexed_inv = DataObject('inv_Jac', *self.tranformation_indices)
        else:
            raise NotImplementedError("")

        remove = lambda x:self.remove_fn(x)
        transformed_vector = zeros(*vector.shape)
        for no,v, in enumerate(vector):
            fd = self.evaluations[order].transformed_derivative[direction]
            # Substitute the DataObject atoms with their corresponding symbols
            fn = v.subs(dict([(s,s.subs(self.cartesian_to_curvilinear_mapping)) for s in v.atoms(DataObject)]))
            # Transform the vector
            fd = fd.subs(DataObject(str(self.wildfunc)), fn)
            transformed_vector[no] = fd
        return transformed_vector

    def transform_first_derivative_scalar(self, scalar, direction):
        order = 1
        if order in self.evaluations.keys():
            if not self.evaluations[order].used:
                self.evaluations[order].used = True
                self.evaluations[order].Jacobian_indexed = DataObject('Jac', *self.tranformation_indices)
                self.evaluations[order].Jacobian_indexed_inv = DataObject('inv_Jac', *self.tranformation_indices)
        else:
            raise NotImplementedError("")

        remove = lambda x:self.remove_fn(x)
        fd = self.evaluations[order].transformed_derivative[direction]
        # Substitute the DataObject atoms with their corresponding symbols
        fn = scalar.subs(dict([(s,s.subs(self.cartesian_to_curvilinear_mapping)) for s in scalar.atoms(DataObject)]))
        # Transform the vector
        fd = fd.subs(CoordinateObject(str(self.wildfunc)), fn)
        return fd

    def transform_second_derivative_scalar(self, scalar):
        """ Transforms the second derivative """
        order == 2
        if order in self.evaluations.keys():
            self.evaluations[order].used = True
        else:
            raise NotImplementedError("")
        return


    def create_transform_equations(self, orders):
        fd_dict = {}
        sd_dict = {}
        fd_eq = []
        sd_eq = []

        if 1 in orders:
            symb = []
            evals = []
            eval_object = self.evaluations[1]
            for i, value in enumerate(eval_object.symbols):
                if value != 0 and isinstance(value,DataObject):
                    symb.append(value)
            for i, value in enumerate(eval_object.metric_evaluations/eval_object.detJ):
                if value !=0:
                    evals.append(value)
            for i in range(len(symb)):
                fd_eq.append(Eq(symb[i], evals[i]))
                fd_dict[symb[i]] = evals[i]
        if 2 in orders:
            symb = []
            evals = []
            eval_object = self.evaluations[2]
            for i, value in enumerate(eval_object.symbols):
                if value != 0:
                    symb.append(value)
            subs_fd = lambda x:x.subs(fd_dict)
            eval_object.eval_symbols = eval_object.eval_symbols.applyfunc(subs_fd)

            for i, value in enumerate(eval_object.eval_symbols):
                if value != 0:
                    evals.append(value)
            for i in range(len(symb)):
                sd_eq.append(Eq(symb[i], evals[i]))
                sd_dict[symb[i]] = evals[i]
        # self.print_latex(sd_eq)
        # pprint(sd_eq)
        return

    def check_einstein(self, equations):
        variables = []
        syms = []
        indexed = []
        for eq in equations:
            variables = variables + list(eq.atoms(CoordinateObject))
            syms = syms + list(eq.atoms(Symbol))
            indexed = indexed + list(eq.atoms(DataObjectBase))
        pprint(set(variables))
        pprint(set(syms))
        pprint(set(indexed))
        pprint(set(variables+indexed).difference(set(syms)))
        pprint([srepr(i) for i in set(indexed)])
        pprint([srepr(i) for i in set(syms)])
        return


    def create_transformation_evaluations(self, spatial_derivative):
        """ Creates the computations Kernel for the evaluation of the coordinate transformation term. This would also populate
        grid.known as well  by default all the derivatives are not stored here they are evaluated as local derivatives
        """
        self.metrics_equiatons = []
        if not bool(self.evaluations):
            raise ValueError("The metric evaluations are not populated with te transformation class\
                intilialized")
        equations = []
        known = []
        self.metric_derivatives_required = set()
        for key, value in self.evaluations.iteritems():
            # Add to the known values as these are computed , this is not true
            # Equation for the determinant of the Jacobian
            if value.used:
                self.metric_derivatives_required = value.detJ.atoms(Derivative).union(self.metric_derivatives_required)
                equations += [Eq(value.Jacobian_indexed, value.detJ)]
                equations += [Eq(value.Jacobian_indexed_inv, 1.0/(value.detJ))]
                known += [d.args[0] for d in value.detJ.atoms(Derivative)] # We assume we know x0...
                if value.order == 1:
                    for number,v in enumerate(list(value.metric_evaluations)):
                        # The equations for the metric terms = adjointJ
                        eq = Eq(value.symbols[number],v/value.detJ)
                        # Add the derivatives in the adjoint Matrix to store as well
                        self.metric_derivatives_required = v.atoms(Derivative).union(self.metric_derivatives_required)
                        known += [d.args[0] for d in v.atoms(Derivative)]
                        #self.metric_derivatives_required
                        if isinstance(eq, Equality):
                            equations +=[eq]
                else:
                    raise NotImplementedError("Not implemented higer derivatives metric eval")
        #pprint(equations)
        val = self.store_derivatives
        self.store_derivatives = False
        from .evaluations import create_formula_evaluations, set_range_of_evaluations, Evaluations
        self.known += list(set(known))
        #pprint(self.known)
        evaluations = {}
        for val in self.known:
            evaluated = Evaluations(val, val, None, None, val)
            evaluations[val] = evaluated
        #evaluations = {}
        self.metrics_equiatons += equations
        #pprint(self.metrics_equiatons)
        pprint([(s, [srepr(s1) for s1 in s.atoms(DataObjectBase)]) for s in self.metrics_equiatons])
        evaluations = spatial_derivative.create_evaluations(self.metrics_equiatons, evaluations, self)
        # TODO update the work arrays for the metric evaluations
        # set the range of evaluations
        #
        #pprint(spatial_derivative.order_of_evaluations)
        set_range_of_evaluations(spatial_derivative.order_of_evaluations, evaluations, self)

        # create the terms
        #evaluations = create_formula_evaluations(self.metrics_equiatons, evaluations)
        #pprint(evaluations)
        #pprint(equations)

        from .utils import update_work_arrays, substitute_work_arrays
        evaluations = update_work_arrays(spatial_derivative.order_of_evaluations,evaluations,self )
        # Create the computations
        computations = []
        # This will return none as we are not saving any of the derivatives
        computations += spatial_derivative.create_computations(self.metrics_equiatons, evaluations, self)
        # Substitute the work arrays in the equations
        evaluation_range = [tuple([0, s]) for s in self.shape]
        updated_equations = substitute_work_arrays(spatial_derivative.order_of_evaluations, evaluations, self.metrics_equiatons, self)
        from .kernel import Kernel
        computations.append(Kernel(updated_equations, evaluation_range, "first_derivative_metric", self))
        # append the LHS of the update equations to Known
        from .grid import GridVariable
        for eq in updated_equations:
            if not isinstance(eq.lhs,GridVariable):
                self.known += [eq.lhs]
        #pprint(self.known)
        # TODO for the boundary conditions we need to apply BC for the first derivative so that it can be used for 2nd
        self.computations = computations
        # similarly do the second derivative metrics
        self.store_derivatives = val
        return

    def transformation_evaluations_BCs(self):
        variables = set()
        for key, value in self.evaluations.iteritems():
            # Add to the known values as these are computed , this is not true
            # Equation for the determinant of the Jacobian
            if value.used:
                variables.add(value.Jacobian_indexed)
                variables.add(value.Jacobian_indexed_inv)
                if value.order == 1:
                    # Order specific evaluated values
                    for number,v in enumerate(list(value.metric_evaluations)):
                        eq = Eq(value.symbols[number],v/value.detJ)
                        if isinstance(eq, Equality):
                            variables.add(value.symbols[number])
                else:
                    raise NotImplementedError("Not implemented higer derivatives metric eval")
        from .bcs import ExtrapolateMetric
        bc = ExtrapolateMetric(self)
        for dim in range(len(self.shape)):
            bc.apply(arrays=variables, boundary_direction=dim ,side=1)
            bc.apply(arrays=variables, boundary_direction=dim ,side=0)
        #(*[comp for comp in bc.computations if comp])
        self.computations += [comp for comp in bc.computations if comp]
        return



class TransformationEvaluation(object):
    """ Object to store the transformation formulae for first and second order transformations."""
    def __init__(self, order):
        self.order = order;
        self.metric_evaluations = None; #What ever the evaluations are (adjA/det(A))..
        self.detJ = None
        self.symbols = None; # The symbols that are to be used
        self.eval_symbols = None
        self.transformed_derivative = None; #The derivative that is transformed
        self.used = False; #This should be True when the transformation is used
        return


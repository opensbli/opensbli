from sympy import Sum, Symbol, Function, flatten, S, pprint, Mul, Expr, Eq
from sympy.tensor import IndexedBase, Indexed, get_contraction_structure
from sympy.tensor.index_methods import get_indices as get_indices_sympy
from sympy import Derivative, preorder_traversal, postorder_traversal
from sympy.tensor.index_methods import _remove_repeated
from opensbli.core.opensbliobjects import DataObject, CoordinateObject, ConstantObject, EinsteinTerm, DataSet, DataSetBase
from sympy.functions.elementary.piecewise import ExprCondPair, Piecewise
from numpy import ndindex as mutidimindex
from sympy.functions.special.tensor_functions import eval_levicivita
from opensbli.utilities.helperfunctions import get_inverse_deltas


class KD(Function):
    """ Handler for the built-in SymPy KroneckerDelta function. """

    @property
    def is_commutative(self):
        return False

    def structure(self):
        """No Contraction structure for KroneckerDelta, returns a indexed object with
        self.args that have a index
        """
        # Create a indexed base
        indexed = IndexedBase("%s" % self.__class__.__name__)
        # Get the arguments to the function that has an index
        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        indexed = indexed[indices]
        indexed.expression = self  # This is for reverse substitution, after einstein expansion
        indexed.is_commutative = False
        return indexed

    def value(self):
        """Evaluates the numerical value of KroneckerDelta"""
        if len(self.args) != 2:
            raise ValueError("Expected only two arguments in KD.")
        if Symbol(str(self.args[0])) == Symbol(str(self.args[1])):
            return S.One
        else:
            return S.Zero


class LC(Function):
    """ Handler for the built-in SymPy LeviCivita function. """

    @property
    def is_commutative(self):
        return False

    def structure(self):
        """No Contraction structure for LeviCivita, returns a indexed object with
        self.args that have a index, see KroneckerDelta above
        """
        if len(self.args) != 3:
            raise ValueError("LeviCivita function should have only three indices.")

        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        indexed_base = IndexedBase("%s" % self.__class__.__name__)
        indexed = indexed_base[tuple(indices)]
        indexed.expression = self
        indexed.is_commutative = False
        return indexed

    def value(self):
        """Evaluates the numerical value of LeviCivita"""
        args = [int(str(a)) for a in self.args]
        return eval_levicivita(*args)


def convert_to_summation(final_expr, sumindices):
    """Convert a given expression into summation notation.
    :param final_expr: the expression to be converted to a SUM.
    :param sumindices: indices over which the expression to be summed
    :returns: Symbolic expression in Summation notation"""
    for k in sumindices:
        final_expr = Sum(final_expr, (k, 0, Symbol('ndim', integer=True)))
    return final_expr


def expand_expr(expr, ndim):
    """Expands a given expression in summation notation to the number of dimensions,
    by applying the Sum

    :param expr: the expression in summation notation to be expanded.
    :param ndim: number of dimensions it should be expanded.
    :returns: expanded expression"""

    def _expand(input_expr, index, ndim):
        """Helper function to expand"""
        expanded_expr = 0
        for dim in range(ndim):
            dictionary = {}
            for at in input_expr.atoms(EinsteinTerm):
                dictionary[at] = at.apply_index(k, dim)
            expanded_expr += input_expr.xreplace(dictionary)
        return expanded_expr

    # Expansion
    # Start with the input expression and apply summation
    expanded_expr = expr.function
    for k in expr.variables:
        expanded_expr = _expand(expanded_expr, k, ndim)
    return expanded_expr


def expand_free_indices(expr, indices, ndim):
    """Expands the free indices, i.e they are not Summed in Einstein expansion.
    For example: In a vector $u_i$, the free index is $i$ this would be expanded
    to $\left[u0, u1, u2\right]$ for ndim 3.

    :param expr: the expression to which free indices are to be applied.
    :param indices: the indices that are to be applied
    :param ndim: number of dimensions it should be expanded.
    :returns: a list of expanded expressions"""

    # iterator for all combinations of the numerical value of indices
    mdindex = mutidimindex(tuple([ndim for i in indices]))  # use numpy's function

    # Create a list to store the output
    # TODO V2 the output can be returned as a list (vector), Matrix for Rank 2 tensor, or sympy
    # multi-dim array for higher order tensors

    expanded_expressions = []
    indices = list(indices)
    for index in mdindex:  # numerical indices combinations
        local_expr = expr.copy()  # copy the given expression
        # create a map of index and value
        local_map = [(indices[no], index[no]) for no, i in enumerate(indices)]
        for at in local_expr.atoms(EinsteinTerm):
            # Apply the mapping to get the component term
            evaluated = at.apply_multiple_indices(indices, local_map)
            local_expr = local_expr.replace(at, evaluated)
        expanded_expressions += [local_expr]
    return expanded_expressions


class EinsteinStructure(object):
    # TODO V2 comments
    def _structure(cls, expr):
        substits = {}
        obj_substitutions = {}
        pot = preorder_traversal(expr)
        # Find the structure of individual local functions in the given expression
        for p in pot:
            if isinstance(p, localfuncs):
                if p.structure():
                    substits[p] = p.structure()
                    expr = expr.xreplace({p: substits[p]})
                else:
                    pass
                pot.skip()
            elif isinstance(p, local_objects):
                if p.structure():
                    obj_substitutions[p] = p.structure()
                else:
                    pass
                pot.skip()
            elif(isinstance(p, Derivative)):
                pprint(p)
                raise ValueError("This derivative is not classified", p)
            else:
                continue
        expr = expr.xreplace(obj_substitutions)
        contraction_dictionary = (get_contraction_structure(expr))
        indices, dummies = get_indices_sympy(expr)
        expr = cls.apply_contraction(contraction_dictionary, expr)
        return expr, indices

    def structure(cls):
        """This finds the Einstein structure for the class by repeatedly calling _structure for each
        argument to the class. _structure returns the summation notation"""
        local_structure = []
        for arg in cls.args:
            expr, outer_indices = cls._structure(arg)
            local_structure += [(expr, outer_indices)]
        all_outer_inds = []
        expr_args = []
        for ind in local_structure:
            a, b = ind
            expr_args += [a]
            all_outer_inds += list(b)
        final_expr = type(cls)(*expr_args)
        if all_outer_inds:
            temporary_indexed_object = IndexedBase("temp")[all_outer_inds]
            inds, summations = _remove_repeated(temporary_indexed_object.indices)
        else:
            summations = None
            inds = None
        if summations:
            final_expr = convert_to_summation(final_expr, summations)
        if inds:
            indexedobj = IndexedBase("%s%s" % (cls.simple_name, cls.hashs(inds)))[tuple(inds)]
        else:
            indexedobj = IndexedBase("%s%s" % (cls.simple_name, cls.hashs()))
        final_expr = cls.substitute_indexed(final_expr)
        indexedobj.expression = final_expr
        return indexedobj

    def hashs(cls, outer_structure=None):
        if outer_structure:
            tup = (cls.args + tuple(list(outer_structure)))
        else:
            tup = cls.args
        return hash((type(cls).__name__,) + tup)

    def apply_contraction(cls, contraction_dictionary, expr):
        outer = []
        replacements = {}
        # Traverse the contraction structure

        def _contraction_traverse(d1, outer, replacements):
            for key in d1:
                if isinstance(key, Expr):
                    continue
                for term in d1[key]:
                    if term in d1:
                        for term1 in d1[term]:
                            _contraction_traverse(term1, outer, replacements)
                if key:
                    for val in d1[key]:
                        summation_notation = val
                        summation_notation = convert_to_summation(summation_notation, key)
                        if val in replacements:
                            raise ValueError("The replacement object already exist this might cause errors")
                        summation_notation = summation_notation.xreplace(replacements)
                        replacements[val] = summation_notation
                    outer += [tuple([key, d1[key]])]

            return replacements, outer
        replacements, outer = _contraction_traverse(contraction_dictionary, outer, replacements)
        expr = expr.subs(replacements)
        return expr

    def substitute_indexed(cls, indexedexpr):
        """Takes an indexed expression and substitutes the original expression for that
        indexed expression"""
        pot = preorder_traversal(indexedexpr)
        indexed_subs = {}
        indexedBaseSubs = {}
        for p in pot:
            if isinstance(p, Indexed):
                indexed_subs[p] = p.expression
                pot.skip()
            elif isinstance(p, IndexedBase):
                indexedBaseSubs[p] = p.expression
                pot.skip()
            else:
                continue
        # now substitute them to get the original expression
        indexedexpr = indexedexpr.xreplace(indexed_subs)
        indexedexpr = indexedexpr.xreplace(indexedBaseSubs)
        return indexedexpr

    def expand_summations(cls, expression, ndim):
        # To expand the summation notation. Do a postorder traversal of the tree
        pot = postorder_traversal(expression)
        to_expand = []
        for p in pot:
            if isinstance(p, Sum):
                to_expand += [p]
        # Once we have the terms to be expanded
        expansion_subs = {}
        for no, val in enumerate(to_expand):
            value = val.subs(expansion_subs)
            if isinstance(value, Sum):
                expansion_subs[val] = expand_expr(value, ndim)
            else:
                expansion_subs[val] = value
        # Do a pre_order traversal and substitute the values in the original expression
        pot = preorder_traversal(expression)
        for p in pot:
            if isinstance(p, Sum):
                expression = expression.xreplace({p: expansion_subs[p]})
                pot.skip()
            else:
                continue
        return expression


class BasicDiscretisation(EinsteinStructure):
    # TODO V2 comments
    @property
    def required_datasets(cls):
        objs = list(cls.args[0].atoms(DataSet))
        return objs

    @property
    def required_datasetbases(cls):
        objs = list(cls.args[0].atoms(DataSetBase)) + list(cls.work.atoms(DataSetBase))
        return list(set(objs))

    @property
    def required_constants(cls):
        constants = list(cls.args[0].atoms(ConstantObject))
        return constants

    @property
    def required_functions(cls):
        """
        # TODO V2 where is it used?
        this should return the functions if any in the expression
        """
        subeval = []
        if cls.args[0].atoms(Function):
            raise ValueError("Argument zero of derivative contains nested ders")
        if cls.is_homogeneous:
            pass
        else:
            subclasstypes = len(cls.args) - 2
            for n in range(subclasstypes):
                newargs = [cls.args[0]] + [cls.args[n+1]]
                subeval += [type(cls)(*newargs)]
        return subeval

    @property
    def order(cls):
        return len(cls.args)-1

    @property
    def is_homogeneous(cls):
        if len(set(cls.args[1:])) == 1:
            return True
        else:
            return False

    @property
    def is_store(cls):
        if hasattr(cls, 'store'):
            return cls.store
        else:
            return False

    @property
    def get_direction(cls):
        direction = []
        for args in cls.args[1:]:
            direction += [args.direction]
        return direction

    @property
    def _sanitise(cls):
        """# This would apply only to derivative types with more than one variable"""
        if not cls.is_homogeneous:
            coordinate_directions = cls.get_direction
            sorted_variables = [x for _, x in sorted(zip(coordinate_directions, cls.args[1:]))]
            args = [cls.args[0], sorted_variables[0]]
            expr = type(cls)(*args)
            for arg in sorted_variables[1:]:
                args = [expr, arg]
                expr = type(cls)(*args)
            return expr
        else:
            return cls

    @property
    def used_else_where(self):
        # TODO V2 : Check if required
        if hasattr(self, '_is_used'):
            return self._is_used
        else:
            return False

    def is_used(self, value):
        # TODO V2 : Check if required
        self._is_used = value
        return

    def update_work(cls, block):
        """ #TODO check if used
        For this we set the indices of the work array to be same as that of
        the derivative indices. This way no mapping is required.
        Block work array will be implemented
        """
        cls.work = block.work_array()
        block.increase_work_index
        return

    def set_store(cls, value):
        cls.store = value
        return


class Dot(Function, BasicDiscretisation):
    # TODO V2 documentation
    def __new__(cls, arg1, arg2):
        ret = super(Dot, cls).__new__(cls, arg1, arg2)
        return ret

    @property
    def simple_name(self):
        return "Dot"

    def value(self):
        return Mul(*self.args)


class DerPrint(object):
    def _pretty(self, printer, *args):
        """Pretty printing class for the derivative"""
        from sympy.printing.pretty.stringpict import prettyForm, stringPict
        from sympy.printing.pretty.pretty_symbology import U
        from sympy.utilities import group

        if printer._use_unicode:
            deriv_symbol = U('PARTIAL DIFFERENTIAL')
        else:
            deriv_symbol = r'd'
        syms = list(reversed(self.args[1:]))
        x = None

        for sym, num in group(syms, multiple=False):
            s = printer._print(sym)
            ds = prettyForm(*s.left(deriv_symbol))

            if num > 1:
                ds = ds**prettyForm(str(num))

            if x is None:
                x = ds
            else:
                x = prettyForm(*x.right(' '))
                x = prettyForm(*x.right(ds))

        f = prettyForm(
            binding=prettyForm.FUNC, *printer._print(self.args[0]).parens())

        pform = prettyForm(deriv_symbol)

        if len(syms) > 1:
            pform = pform**prettyForm(str(len(syms)))

        pform = prettyForm(*pform.below(stringPict.LINE, x))
        pform.baseline = pform.baseline + 1
        pform = prettyForm(*stringPict.next(pform, f))
        pform.binding = prettyForm.MUL
        return pform

    def _sympystr(self, p):
        """This is the string used as a representation without any `commas`"""
        args = list(map(p.doprint, self.args))
        return "%s %s" % (self.simple_name, " ".join(args))


class CentralDerivative(Function, BasicDiscretisation, DerPrint):
    """
    TODO V2 documentation
    wrapper class to represent derivatives
    Sympy already have a "Derivative" class, thus double D
    """
    def __new__(cls, expr, *args):
        args = tuple(flatten([expr] + list(args)))
        ret = super(CentralDerivative, cls).__new__(cls, *args, evaluate=False)
        ret.store = True  # By default all the derivatives are stored
        ret.local_evaluation = True
        return ret

    def doit(cls):
        if any(arg == S.Zero for arg in cls.args):
            return S.Zero
        elif len(set(cls.args)) == 1:
            return S.One
        else:
            return cls

    def expand(self, **hints):
        from sympy.core.function import _coeff_isneg
        ders = self.args[1:]
        rets = 0
        arg = self.args[0]
        if arg.is_Add:
            aargs = list(arg.expand(deep=True).args)
            for i, ai in enumerate(aargs):
                if ai.is_Mul and _coeff_isneg(ai):
                    ai = ai*-S.One
                    rets -= self.func(ai, *ders)
                else:
                    rets += self.func(ai, *ders)
            return rets
        else:
            return self

    def _eval_expand_func(self, **hints):
        return self.expand()

    def _discretise_derivative(cls, scheme, block, boundary=True):
        """
        TODO V2 documentation
        This would return the discritized derivative of the
        local object depending on the order of accuracy specified
        Returns the formula for the derivative function, only first derivatives or homogeneous
        derivatives of higher order are supported. The mixed derivatives will be handled impl-
        citly while creating the kernels
        :arg derivative: the derivative on which discretisation should be performed
        :returns: the discritized derivative, in case of wall boundaries this is a Piecewise-
        function
        """
        order = cls.order
        form = 0
        # Put the coefficients of first and second derivatives in a dictionary and use them
        if cls.is_homogeneous:
            dire = cls.get_direction[0]
            weights = scheme._generate_weights(dire, order, block)
            for no, p in enumerate(scheme.points):
                expr = cls.args[0]
                for req in (cls.required_datasets):
                    loc = list(req.indices)
                    loc[dire] = loc[dire] + p
                    val = req.base[loc]
                    expr = expr.replace(req, val)
                form = form + weights[no]*expr
            if form == 0:
                raise ValueError("Central derivative formula is zero for %s" % cls)
        else:
            raise ValueError("The provided derivative is not homogeneous, %s" % cls)
        if boundary:
            form = cls.modify_boundary_formula(form, block)

        delta = S.One/block.deltas[dire]**order
        inv_delta = get_inverse_deltas(delta)
        form = form*(inv_delta)
        return form

    def modify_boundary_formula(cls, form, block):
        # Apply the boundary modifications
        modifications = block.check_modify_central()
        dire = cls.get_direction[0]
        if dire in modifications:
            boundary_mods = [k for k in modifications[dire] if k]
            expression_condition_pairs = []
            for b in boundary_mods:
                expression_condition_pairs += b.modification_scheme.expr_cond_pairs(cls.args[0], b.direction, b.side, cls.order, block)
            expression_condition_pairs += [ExprCondPair(form, True)]
            form = Piecewise(*expression_condition_pairs, **{'evaluate': False})
        return form

    @property
    def simple_name(cls):
        return "%s" % ("CD")

    def classical_strong_differentiabilty_transformation(cls, metric):
        direction = cls.get_direction
        if cls.order == 1:
            metric_der = metric.classical_strong_differentiabilty_transformation[direction[0]]
            transformed_der = metric_der.subs(metric.general_function, cls.args[0])
        elif cls.order == 2:
            metric.sd_used = True
            metric_der = metric.classical_strong_differentiabilty_transformation_sd[tuple(cls.get_direction)]
            transformed_der = metric_der.subs(metric.general_function, cls.args[0])
        return transformed_der


class WenoDerivative(Function, BasicDiscretisation, DerPrint):

    def __new__(cls, expr, *args, **settings):
        args = flatten([expr] + list(args))
        ret = super(WenoDerivative, cls).__new__(cls, *args, evaluate=False)
        ret.store = True  # By default all the derivatives are stored
        ret.reconstructions = []
        ret.local_evaluation = True
        ret.settings = settings
        return ret

    @property
    def simple_name(cls):
        return "%s" % ("WD")

    def update_settings(self, **settings):
        # existing_keys = self.settings.keys()
        for key in settings.keys():
            if key in self.settings.keys():
                raise ValueError("Key exists")
        self.settings.update(settings)
        return

    def _discretise_derivative(cls, block, scheme=None):
        """This would return the discritized derivative of the
        local object depending on the order of accuracy specified
        Returns the formula for the derivative function, only first derivatives or homogeneous
        derivatives of higher order are supported. The mixed derivatives will be handled impl-
        citly while creating the kernels
        :arg derivative: the derivative on which discretisation should be performed
        :returns: the discritized derivative, in case of wall boundaries this is a Piecewise-
        function
        """
        order = cls.order
        if (order > 1):
            raise ValueError("Weno Derivatives only defined for first order")
        dire = cls.get_direction[0]
        delta = block.deltas[dire]
        loc = list(cls.reconstruction_work.indices[:])
        loc[dire] += -1
        form = (cls.reconstruction_work - cls.reconstruction_work.base[loc]) / delta
        return form

    def add_reconstruction_classes(self, classes):
        if isinstance(classes, list):
            self.reconstructions += classes
        else:
            self.reconstructions += [classes]
        return

    def create_reconstruction_work_array(self, block):
        self.reconstruction_work = block.work_array()
        block.increase_work_index
        return

    def classical_strong_differentiabilty_transformation(cls, metric):
        direction = cls.get_direction
        if cls.order == 1:
            metric_der = metric.classical_strong_differentiabilty_transformation[direction[0]]
        elif cls.order == 2:
            raise NotImplementedError("")
        for at in metric_der.atoms(Function):
            local_at = type(cls)(at.args[0].subs(metric.general_function, cls.args[0]), at.args[1:])
            metric_der = metric_der.subs(at, local_at)
        return metric_der

    @property
    def evaluate_reconstruction(self):
        # Check if the reconstruction placeholders are combined, other options should be added here
        if "combine_reconstructions" in self.settings and self.settings["combine_reconstructions"]:
            variables = set([r.reconstructed_symbol for r in self.reconstructions])
            if len(variables) == 1:
                return list(variables)[0]
            else:
                raise ValueError("")
        else:
            # Sum the reconstruction placeholders
            total = 0
            for r in self.reconstructions:
                total += r.reconstructed_symbol
            return total


class TenoDerivative(Function, BasicDiscretisation, DerPrint):

    def __new__(cls, expr, *args, **settings):
        args = flatten([expr] + list(args))
        ret = super(TenoDerivative, cls).__new__(cls, *args, evaluate=False)
        ret.store = True  # By default all the derivatives are stored
        ret.reconstructions = []
        ret.local_evaluation = True
        ret.settings = settings
        return ret

    @property
    def simple_name(cls):
        return "%s" % ("TD")

    def update_settings(self, **settings):
        # existing_keys = self.settings.keys()
        for key in settings.keys():
            if key in self.settings.keys():
                raise ValueError("Key exists")
        self.settings.update(settings)
        return

    def _discretise_derivative(cls, block, scheme=None):
        """This would return the discritized derivative of the
        local object depending on the order of accuracy specified
        Returns the formula for the derivative function, only first derivatives or homogeneous
        derivatives of higher order are supported. The mixed derivatives will be handled impl-
        citly while creating the kernels
        :arg derivative: the derivative on which discretisation should be performed
        :returns: the discritized derivative, in case of wall boundaries this is a Piecewise-
        function
        """
        order = cls.order
        if (order > 1):
            raise ValueError("Teno Derivatives only defined for first order")
        dire = cls.get_direction[0]
        delta = block.deltas[dire]
        loc = list(cls.reconstruction_work.indices[:])
        loc[dire] += -1
        form = (cls.reconstruction_work - cls.reconstruction_work.base[loc]) / delta
        return form

    def add_reconstruction_classes(self, classes):
        if isinstance(classes, list):
            self.reconstructions += classes
        else:
            self.reconstructions += [classes]
        return

    def create_reconstruction_work_array(self, block):
        self.reconstruction_work = block.work_array()
        block.increase_work_index
        return

    @property
    def evaluate_reconstruction(self):
        # Check if the reconstruction placeholders are combined, other options should be added here
        if "combine_reconstructions" in self.settings and self.settings["combine_reconstructions"]:
            variables = set([r.reconstructed_symbol for r in self.reconstructions])
            if len(variables) == 1:
                return list(variables)[0]
            else:
                raise ValueError("")
        else:
            # Sum the reconstruction placeholders
            total = 0
            for r in self.reconstructions:
                total += r.reconstructed_symbol
            return total

    def classical_strong_differentiabilty_transformation(cls, metric):
        direction = cls.get_direction
        if cls.order == 1:
            metric_der = metric.classical_strong_differentiabilty_transformation[direction[0]]
        elif cls.order == 2:
            raise NotImplementedError("")
        for at in metric_der.atoms(Function):
            local_at = type(cls)(at.args[0].subs(metric.general_function, cls.args[0]), at.args[1:])
            metric_der = metric_der.subs(at, local_at)
        return metric_der


class TemporalDerivative(Function, BasicDiscretisation, DerPrint):

    def __new__(cls, expr, *args):
        args = flatten([expr] + list(args))
        ret = super(TemporalDerivative, cls).__new__(cls, *args, evaluate=False)
        ret.store = True  # By default all the derivatives are stored
        ret.local_evaluation = True
        return ret

    @property
    def simple_name(cls):
        return "%s" % ("TimeDer")

    @property
    def time_advance_array(cls):
        return cls.args[0]

    def classical_strong_differentiabilty_transformation(cls, metric):
        # No change to the class
        return cls


class MetricDerivative(Function, BasicDiscretisation):
    """
    # TODO V2: Documentation
    """
    def __new__(cls, expr, *args):
        args = tuple(flatten([expr] + list(args)))
        ret = super(MetricDerivative, cls).__new__(cls, *args, evaluate=False)
        ret.store = True  # By default all the derivatives are stored
        return ret


localfuncs = (MetricDerivative, KD, CentralDerivative, WenoDerivative, TenoDerivative, TemporalDerivative, LC, Dot)
simplifying_funcs = (KD, LC, Dot)
local_objects = (DataObject, CoordinateObject, ConstantObject, EinsteinTerm)

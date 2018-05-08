from sympy import Symbol, flatten
from sympy.core.compatibility import is_sequence
from sympy.tensor import Idx, IndexedBase, Indexed
from sympy import pprint, srepr
from sympy.tensor.indexed import IndexException
from sympy.core.cache import cacheit
from opensbli.core.datatypes import SimulationDataType
from sympy.core import Basic, Tuple, Function, Equality
from sympy.core.basic import *
from sympy.logic.boolalg import Boolean


_projectname = "opensbli"

class EinsteinTerm(Symbol):

    """ Basic term in OpenSBLI. This is used to parse the equations. Any atom in the expression that is not a function
    is represented as an EinsteinTerm. Atoms like constant, coordinates, Data objects are derived from this.

    .. note::
        This could be e.g. :math:`{\\tau}_{ij}`, but can also be e.g. :math:`{u}_{i}, \\rho, {x}_{j}` (coordinates), Re (constants).
        All symbols in the equation are Einstein terms, but they can have zero or more indices. The indices
        are used for applying the Einstein contraction structure.
        **Used during parsing and Einstein expansion process**

    :arg str name: The name of the symbol under consideration. This can have zero or more indices."""

    is_commutative = False
    
    def __new__(self, name, **assumptions):
        self._sanitize(assumptions, self)  # Remove any 'None's, etc.
        self.name = str(name)

        # Make this into a new SymPy Symbol object.
        self = Symbol.__xnew__(self, self.name, **assumptions)

        # Is this a term that is constant in space and time (i.e. doesn't have any indices).
        self.is_constant = False
        self.is_coordinate = False
        self.is_tensor = False
        self.rank = 0  # Rank of the term

        # Extract the indices, which are always preceded by an underscore.
        indices = self.name.split('_')[1:]
        self.indices = [Idx(x) for x in indices]
        if self.indices:
            self.is_tensor = True
            self.rank = len(self.indices)
        return self

    def get_indices(self):
        """Return a list of the Einstein indices.

        :returns: A list of the Einstein indices.
        :rtype: list """
        return self.indices

    def structure(self):
        """Structure of the Einstein term, used for determining the contraction structure of an expression

        :returns: SymPy's Indexed object with Einstein term base as the base and indices are Einstein term indices
        :rtype: Indexed """
        if len(self.get_indices()) > 0:
            st = IndexedBase("%s" % str(self.get_base()))
            st = st[self.get_indices()]
            st.expression = self
            self.indexed_object = st
            return st
        else:
            st = self
            st.expression = st
            self.indexed_object = st
            return st

    def get_base(self):
        """Return the base name of the EinsteinTerm without indices.

        :returns: The base name.
        :rtype: str """
        return self.name.split('_')[0]

    def apply_index(self, idx, value):
        """Substitutes the index *idx* with its value and returns a new Einstein term.

        :param Idx idx: Index of the original term to be substituted
        :param int value: value of the index
        :returns: a new Einstein term with the value substituted
        :rtype: EinsteinTerm """
        if idx in self.get_indices():
            if self.get_base():
                newval = str(self).replace("_%s" % str(idx), str(value))
                val = type(self)(newval)
                indices = list(set(self.get_indices()[:]).difference(set([idx])))
                val.indices = indices
                # Store val index
                val.direction = value
            else:
                val = int(value)
        else: # Idx doesn't exist in the EinsteinTerm.
            val = self
        return val

    def apply_multiple_indices(self, indices, index_map):
        """Similar to apply_index function but can handle multiple indices.

        :param indices: Indices to be applied.
        :param index_map: A map of the index and its value.
        :returns: A new Einstein term with all of the index values substituted.
        :rtype: EinsteinTerm
        """
        dictionary_indices = {}
        newval = str(self)
        for m in index_map:
            newval = newval.replace("_%s" % str(m[0]), str(m[1]))
            dictionary_indices[m[0]] = m[1]
        val = type(self)(newval)
        # Store value of the index if it is a CoordinateObject
        if isinstance(self, CoordinateObject):
            if self.get_indices():
                val.direction = dictionary_indices[self.get_indices()[0]]
        return val


class Constant(object):
    """ A base class to represents different types of constants in OpenSBLI."""
    pass


class ConstantObject(EinsteinTerm, Constant):
    """A constant object which can have Einstein indices to be expanded. This is used to
    differentiate between different Einstein terms, which are used in differentiation.

    **Used during parsing and Einstein expansion process**

    :param str label: name of the constant object
    :returns: declared constant
    :rtype: ConstantObject """
    is_commutative = True

    def __new__(cls, label, **kwargs):
        ret = super(ConstantObject, cls).__new__(cls, label, **kwargs)
        ret.is_constant = True
        ret.is_input = True
        ret._datatype = SimulationDataType()
        ret._value = "Input"
        return ret
    
    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return str(self.name)

    @property
    def datatype(self):
        """Numeric datatype of the ConstantObject.

        :returns: Numerical datatype (see :class:`.SimulationDataType`)
        :rtype: str """
        return self._datatype

    @datatype.setter
    def datatype(self, dtype):
        """Set the data type of the Constant."""
        self._datatype = dtype

    @property
    def value(self):
        """Returns the value of Constant."""
        return self._value
    
    @value.setter
    def value(self, numerical_value, dtype=None):
        """Sets the value of the Constant.
        :param numerical_value: Value to be set for the Constant.
        :param dtype: Data type of the Constant, defaults to SimulationDataType. """
        self.is_input = False
        self._value = numerical_value
        if dtype:
            self.datatype = dtype
        else:
            self.datatype = SimulationDataType()
        return


class ConstantIndexed(Indexed, Constant):
    """ An indexed object represented by an array of constants.

    :param str label: Name of the ConstantIndexed.
    :param list indices: Indices of the ConstantIndexed. (See: Sympy Indexed class)."""
    def __new__(cls, label, indices, **kwargs):
        base = IndexedBase(label)
        if isinstance(indices, list):
            for i in indices:
                if not isinstance(i, Idx):
                    raise ValueError("The indices of the Constant Indexed Object should be of type Idx", i)
        else:
            if not isinstance(indices, Idx):
                raise ValueError("The indices of the Constant Indexed Object should be of type Idx", i)
            indices = flatten([indices])
        ret = Indexed.__new__(cls, base, *indices)
        ret.is_constant = True
        ret.inline_array = True
        ret.is_input = True
        ret._datatype = SimulationDataType()
        ret._value = ["Input" for i in range(ret.shape[0])]
        return ret

    @property
    def datatype(self):
        """Numeric data type of the constant array.

        :returns: Numerical datatype (see :class:`.SimulationDataType`)"""
        return self._datatype

    @datatype.setter
    def datatype(self, dtype):
        self._datatype = dtype

    @property
    def value(self):
        return self._value

    @property
    def name(self):
        return str(self.base)
    
    @value.setter
    def value(self, numerical_values, dtype=None):
        self.is_input = False
        if len(numerical_values) != len(self.value):
            raise ValueError("Values for ConstantIndexed should be of length of the constants.")
        self._value = numerical_values
        if dtype:
            self.datatype = dtype
        else:
            self.datatype = SimulationDataType()
        return
    
    @property
    def value_access_c(self):
        """Returns the C code for accessing the values of the ConstantIndexedObject."""
        return ["%s[%d]" % (self.base, i) for i in range(0, len(self.value))]

    @property ## WARNING: Is this required?
    def location(cls):
        return list(cls.args[1:])


class CoordinateObject(EinsteinTerm):
    """ A coordinate object which can have Einstein indices to be expanded, this is used to
    differentiate between different Einstein terms, while performing differentiation.

    **Used during parsing, Einstein expansion processes and during discretisation**

    :param str label: Name of the coordinate object to be defined
    :returns: Declared CoordinateObject.
    :rtype: CoordinateObject """
    is_commutative = True
    is_Atom = True

    @cacheit
    def __new__(cls, label, **kwargs):
        ret = super(CoordinateObject, cls).__new__(cls, label)
        if 'time' in kwargs:
            ret.timecoordinate = True
        else:
            ret.timecoordinate = False
        ret.is_coordinate = True
        return ret

    def get_coordinate_type(cls):
        """Returns True if the CoordinateObject is a time coordinate."""
        if cls.timecoordinate:
            return True
        else:
            return False

    def set_direction(cls, dire):
        cls.direction = dire
        return


class MetricObject(EinsteinTerm):
    """ A metric object which can have Einstein indices to be expanded, this is used to
    differentiate between different Einstein terms, while performing metric transformations

    **Used during application of Metrics and expansion**

    :param str label: name of the metric object to be defined
    :returns: declared metric
    :rtype: MetricObject """
    def __new__(cls, label, **kwargs):
        ret = super(MetricObject, cls).__new__(cls, label)
        if 'time' in kwargs:
            ret.timecoordinate = True
        else:
            ret.timecoordinate = False
        return ret


class DataObject(EinsteinTerm):
    """ Once parsed and expanded all the EinsteinTerms that are not constants and coordinates are converted
    to DataObject. This is useful to discretise the expanded equations on a block and acts as a link between
    Equation expansion and discretisation. If required the user can write their own equations using these objects.

    This acts as an intermediate function to decouple the expansion process and other opensbli processes, and
    gives flexibility in equation definition.

    :param str label: Name of the DataObject to be defined.
    :returns: Declared DataObject.
    :rtype: DataObject """
    is_commutative = True
    is_Atom = True

    def __new__(cls, label, **kw_args):
        ret = super(DataObject, cls).__new__(cls, label, **kw_args)
        return ret

    def copy(self):
        return self

    @property
    def free_symbols(self):
        return {self}


class DataSetBase(IndexedBase):
    """ Base object for converting a DataObject to an array on the block with specific dimensions of
    the block. When a SimulationBlock is instantiated then the attribute block is set and the various
    attributes of the block are used.

    :param label: Name of the Dataset base.
    :returns: The DataSetBase.
    :rtype: DataSetBase

    .. note::
        As given in SymPy's Indexed object documentation it is always advised to create a DataSetBase and then
        a Dataset """
    is_commutative = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True

    @cacheit  # we cache it so that the changes to a particular dataset base are global and not local
    def __new__(cls, label, shape, blocknumber, **kw_args):
        sym = label
        if shape is None:
            raise ValueError("Dataset base requires shape of the block")
        ret = super(DataSetBase, cls).__new__(cls, sym, shape, **kw_args)
        ret.blocknumber = blocknumber
        ret._args = tuple(list(ret._args) + [Idx(blocknumber)])
        return ret

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h
    
    @property
    def blockname(self):
        return "%sblock%02d" %(_projectname, self.blocknumber)

    @property
    def noblockname(self):
        return EinsteinTerm(str(self.label))

    def _hashable_content(self):
        return str(self.label) + self.blockname

    def check_index(cls, indices):
        if len(indices) == len(cls.shape):
            pass
        elif len(indices) != len(cls.shape):
            raise IndexException("Rank mismatch.")
        return indices

    def __getitem__(cls, indices, **kw_args):
        indices = cls.check_index(indices)
        return DataSet(cls, *indices)

    def _sympystr(self, p):
        """ For clarity the block number is printed"""
        return "%s_B%s" % (str(self.label), str(self.blocknumber))

    def simplelabel(self):
        """Returns the base label in case we need the label as a string to modify and create a new one with modifications."""
        return "%s" % (str(self.label))

    @property
    def location(self):
        """ The location is the relative grid location, it is presently hard coded to the grid location.
        This is provided so that if in future staggered grid arrangement can be implemented."""
        return [0 for i in range(len(self.shape))]


class DataSet(Indexed):
    """ It is the Data set of a DatasetBase which is defined on a block.

    :param base: name of the Dataset base, this is of type DataSetBase
    :param list indices: indices of the dataset, in out context it is the value of the dataset at the relative location from the dataset base location (example, [0,0] or [1,0])
    :returns: Dataset at the specified location
    :rtype: DataSet

    .. note::
        As given in SymPy's Indexed object documentation it is always advised to create a DataSetBase and then a Dataset. """
    is_commutative = True
    is_Indexed = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *indices, **kwargs):
        if not isinstance(base, DataSetBase):
            raise ValueError("Declare DatasetBase and instantiate a dataset")
        indices = base.check_index(indices)
        ret = Indexed.__new__(cls, base, *indices)
        return ret

    def _sympystr(self, p):
        return "%s" % (p.doprint(self.base))

    def _pretty(self, printer):
        """Pretty Printing method. """
        from sympy.printing.pretty.stringpict import prettyForm
        allinds = [i for i in self.indices if not isinstance(i, Idx)]
        pforms = [printer._print(a) for a in allinds]
        expr = prettyForm(*printer.join(",", pforms).parens(left="[", right="]"))
        expr = prettyForm(*expr.left(printer._print(self.base)))
        return expr

    @property
    def get_grid_indices(self):
        """ Returns the relative location of the dataset, as used in discretisation
        :rtype: list """
        return [i for i in self.indices if not isinstance(i, Idx)]


class GridIndex(Indexed):
    """ An Indexed object to get the local grid point. """
    is_commutative = True
    is_Indexed = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *indices, **kwargs):
        if not isinstance(base, GridIndexedBase):
            raise ValueError("GridIndex base should be GridIndexedBase object")
        ret = Indexed.__new__(cls, base, *indices)
        return ret
    # def free_symbols(self):
        # return self


class GridIndexedBase(IndexedBase):
    """ Base object to locate the global index of the grid point NOT USED ANY MORE WARNING
    """
    is_commutative = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, label, ndim, **kw_args):
        sym = label
        pprint(sym)
        print type(sym)
        ret = super(GridIndexedBase, cls).__new__(cls, sym, (1), **kw_args)  # Shape would be of size ndim
        pprint(ret.shape)
        return ret

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return str(self.label)

    def __getitem__(cls, indices, **kw_args):
        if isinstance(indices, int):
            indices = [indices]
        if len(indices) != cls.shape:
            raise IndexException("Rank mismatch.")
        return GridIndex(cls, *indices)

    def _sympystr(self, p):
        """ For clarity the block number is printed"""
        return "%s" % (str(self.label))


class Grididx(Symbol):
    """ A coordinate object which can have Einstein indices to be expanded, this is used to
    differentiate between different Einstein terms, while performing differentiation.

    **Used during parsing, Einstein expansion processes and during discretisation**

    :param str label: name of the coordinate object to be defined
    :returns: declared coordinate
    :rtype: CoordinateObject """
    def __new__(self, label, number, **assumptions):
        self._sanitize(assumptions, self)  # Remove any 'None's, etc.
        self.name = str(label) + str(number)
        # Make this into a new SymPy Symbol object.
        self = Symbol.__xnew__(self, self.name, **assumptions)
        self.number = number
        self.base = label
        return self

class GlobalValue(object):
    """A base class for values that are global to all processes, these are not
    constants but they change with time or kind of mesh parameters which are 
    updated in the simulation and not user input"""
    

class Globalvariable(EinsteinTerm, GlobalValue):
    
    is_commutative = True
    
    def __new__(cls, label, **kwargs):
        ret = super(Globalvariable, cls).__new__(cls, label, **kwargs)
        ret._datatype = SimulationDataType()
        ret.is_input = True
        ret._value = "Input"
        return ret
    
    @property
    def datatype(self):
        """Numeric data type of the Globalvariable array.

        :returns: Numerical datatype (see :class:`.SimulationDataType`)"""
        return self._datatype

    @datatype.setter
    def datatype(self, dtype):
        self._datatype = dtype

    @property
    def value(self):
        return self._value


class GroupedCondition(Tuple):
    """Represents an expression, condition pair."""

    def __new__(cls, expressions, condition):
        return Tuple.__new__(cls, expressions, condition)

    @property
    def expressions(self):
        """
        Returns the expression of this pair.
        """
        return self.args[0]

    @property
    def cond(self):
        """
        Returns the condition of this pair.
        """
        return self.args[1]

    @property
    def lhs(self):
        return [eq.lhs for eq in self.expressions]

    @property
    def rhs(self):
        return [eq.rhs for eq in self.expressions]

    @property
    def is_commutative(self):
        return self.expressions.is_commutative

    def __iter__(self):
        yield self.expressions
        yield self.cond
from sympy import Piecewise
class GroupedPiecewise(Piecewise):
    nargs = None
    is_Piecewise = True

    def __new__(cls, *args, **options):
        return Piecewise.__new__(cls, *args, **options)

    def add_pair(cls, expr_pair):
        cls.pairs.append(expr_pair)
        cls.extract_all_equations
        cls.extract_all_conditions
        # assert isinstance(Boolean, expr_pair[1])
        cls.grouped_conditions += [expr_pair[1]]
        cls.grouped_equations += [expr_pair[0]] ## add input checking here
        return
    
    #@property
    #def gro
    
    def convert_to_datasets(self, block):
        replacements = {}
        for d in self.atoms(DataObject):
            replacements[d] = block.location_dataset(d)
        #for expr, c in self.args:
            #if is_sequence(expr):
                #for eq1 in e:
                    #if is_sequence(eq1):
                        #raise NotImplementedError("")
                    #eq1.convert_to_datasets(block)
                
        #for index, list_of_eqn in enumerate(self.grouped_equations):
            #for eqn_no, equation in enumerate(list_of_eqn):
                #self.grouped_equations[index][eqn_no] = equation.convert_to_datasets(block)
        #for index, condition in enumerate(self.grouped_conditions):
            #for d in condition.atoms(DataObject):
                #new = block.location_dataset(str(d))
                #condition = condition.subs(d, new)
            #self.grouped_conditions[index] = condition
        return self.subs(replacements)
    
    def _eval_subs(self, old, new):
        args = list(self.args)
        for i, (e, c) in enumerate(args):
            c = c._subs(old, new)
            if is_sequence(e):
                for no, eq1 in enumerate(e):
                    if is_sequence(eq1):
                        raise NotImplementedError("")
                    e[no] = eq1._subs(old, new)
            else:
                e = e._subs(old, new)
            args[i] = (e, c)
        return self.func(*args)
    
    #@property
    #def required_datasets(self):
        #dsets = set()
        #for eq in flatten(self.grouped_equations):
            #dsets = dsets.union(eq.required_datasets)
        #for c in flatten(self.grouped_conditions):
            #dsets = dsets.union(d.atoms(DataSet))
        #return dsets
    
    @property
    def lhs_datasetbases(self):
        """ These are the datsets to be written out, so only expressions are considered and condition is omitted"""
        dsets = set()
        for e, c in self.args:
            if is_sequence(e):
                for eq1 in e:
                    if is_sequence(eq1):
                        raise NotImplementedError("")
                    dsets = dsets.union(eq1.lhs_datasetbases)
            else:
                dsets = dsets.union(e.lhs_datasetbases)
        return dsets
    
    @property
    def lhs_datasets_full(self):
        """ These are the datsets to be written out, so only expressions are considered and condition is omitted"""
        dsets = set()
        for e, c in self.args:
            if is_sequence(e):
                for eq1 in e:
                    if is_sequence(eq1):
                        raise NotImplementedError("")
                    dsets = dsets.union(eq1.atoms(DataSet))
            else:
                dsets = dsets.union(e.atoms(DataSet))
        return dsets
    
    @property
    def rhs_datasetbases(self):
        dsets = set()
        for e, c in self.args:
            if is_sequence(e):
                for eq1 in e:
                    if is_sequence(eq1):
                        raise NotImplementedError("")
                    dsets = dsets.union(eq1.rhs_datasetbases)
            else:
                dsets = dsets.union(e.rhs_datasetbases)
            dsets = dsets.union(c.atoms(DataSetBase))
        return dsets

    @property
    def expr_rhs(cls):
        return flatten([pairs.rhs for pairs in cls.pairs])

    @property
    def expr_lhs(cls):
        return flatten([pairs.lhs for pairs in cls.pairs])

    @property
    def extract_all_equations(cls):
        cls.all_equations = flatten([pair.expressions for pair in cls.pairs])

    @property
    def extract_all_conditions(cls):
        cls.all_conditions = flatten([pair.cond for pair in cls.pairs])

    
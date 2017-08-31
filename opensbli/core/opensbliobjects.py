from sympy import Symbol, Eq, flatten, srepr
from sympy.tensor import Idx, IndexedBase, Indexed
from sympy import pprint


class EinsteinTerm(Symbol):

    """ Basic term in OpenSBLI. This is used to parse the equations. Any atom in the expression that is not a function
    is represented as an EinsteinTerm. Atoms like constant, coordinates, Data objects are derived from this.

    This could be e.g. tau_i_j, but can also be e.g. u_i, rho, x_j (coordinates), Re (constants).
    In other words, all symbols in the equation are Einstein terms, but they can have zero or more indices. The indices
    are used for applying the Einstein contraction structure.

    **Used during parsing and Einstein expansion process**

    By default EinsteinTerm is not commutative.
    """

    is_commutative = False

    def __new__(self, symbol, **assumptions):
        """ Create a new EinsteinTerm.

        :arg str symbol: The symbol under consideration. This can have zero or more indices.
        """

        self._sanitize(assumptions, self)  # Remove any 'None's, etc.
        self.name = str(symbol)

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
        """ Return a list of the Einstein indices.

        :returns: A list of the Einstein indices.
        :rtype: list
        """
        return self.indices

    def structure(self):
        """ Structure of the Einstein term, used for determining the contraction structure of an expression

        :returns: SymPy's Indexed object with Einstein term base as the base and indices are Einstein term indices
        """
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
        """ Return the base name.

        :returns: The base name.
        :rtype: str
        """
        return self.name.split('_')[0]

    def apply_index(self, idx, value):
        """ Substitutes the index idx with its value and returns a new Einstein term

        :param Idx idx: Index of the original term to be substituted
        :param int value: value of the index
        :returns: a new Einstein term with the value substituted
        :rtype: EinsteinTerm
        """
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
        else:
            val = self
        return val

    def apply_multiple_indices(self, indices, index_map):
        """ Similar to apply_index function but can handle multiple indices

        :param indices: indices to map
        :param index_map: a map of index and its value
        :returns: a new Einstein term with the values substituted
        :rtype: EinsteinTerm
        """
        dictionary_indices = {}
        newval = str(self)
        for m in index_map:
            newval = newval.replace("_%s" % str(m[0]), str(m[1]))
            dictionary_indices[m[0]] = m[1]
        val = type(self)(newval)
        # Store val index
        if isinstance(self, CoordinateObject):
            if self.get_indices():
                val.direction = dictionary_indices[self.get_indices()[0]]
        return val


class Constant(object):
    """ A base class to represents different types of constants
    """
    pass


class ConstantObject(EinsteinTerm, Constant):
    """ A constant object which can have Einstein indices to be expanded. This is used to
    differentiate between different Einstein terms, which are used in differentiation.

    **Used during parsing and Einstein expansion process**

    :param str label: name of the constant object
    :returns: declared constant
    :rtype: ConstantObject
    """
    is_commutative = True

    def __new__(cls, label, **kwargs):
        ret = super(ConstantObject, cls).__new__(cls, label, **kwargs)
        ret.is_constant = True
        return ret


class ConstantIndexed(Indexed, Constant):
    """ A constant Indexed object
    """
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
        return ret

    @property
    def location(cls):
        return list(cls.args[1:])


class CoordinateObject(EinsteinTerm):
    """ A coordinate object which can have Einstein indices to be expanded, this is used to
    differentiate between different Einstein terms, while performing differentiation.

    **Used during parsing, Einstein expansion processes and during discretisation**

    :param str label: name of the coordinate object to be defined
    :returns: declared coordinate
    :rtype: CoordinateObject
    """
    is_commutative = True
    is_Atom = True

    def __new__(cls, label, **kwargs):
        ret = super(CoordinateObject, cls).__new__(cls, label)
        if 'time' in kwargs:
            ret.timecoordinate = True
        else:
            ret.timecoordinate = False
        return ret

    def get_coordinate_type(cls):
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
    :rtype: MetricObject

    """
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
    gives flexibility in equation definition

    :param str label: name of the data objectto be defined
    :returns: declared data object
    :rtype: DataObject

    """
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


from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable, range
from sympy.core.cache import cacheit


class DataSetBase(IndexedBase):
    """ Base object for converting a DataObject to an array on the block with specific dimensions of
    the block. When a SimulationBlock is instantiated then the attribute block is set and the various
    attributes of the block are used.

    :param label: name of the Dataset base, this can be OpenSBLI/SymPy's basic types
    :returns: defined dataset base on the block
    :rtype: DataSetBase which is derived from IndexedBase

    .. note::
        As given in SymPy's Indexed object documentation it is always advised to create a DataSetBase and then
        a Dataset

    """
    is_commutative = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True
    block = None

    @cacheit  # we cache it so that the changes to a particular dataset base are global and not local
    def __new__(cls, label, shape=None, **kw_args):
        if not cls.block:
            raise ValueError("Set the block for DataSetBase")
        sym = label
        if shape is None:
            # print cls.block, "yes"
            shape = list(cls.block.shape) + [Idx(cls.block.blockname)]
        else:
            shape = shape
        ret = super(DataSetBase, cls).__new__(cls, sym, shape, **kw_args)
        ret.noblockname = Symbol(str(label))
        ret.blockname = cls.block.blockname
        ret.blocknumber = cls.block.blocknumber
        return ret

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return str(self.label) + self.blockname

    def __getitem__(cls, indices, **kw_args):
        if len(indices) == len(cls.shape) and indices[-1] == Idx(cls.blockname):
            pass
        elif len(indices) == len(cls.shape) - 1:
            indices += [Idx(cls.blockname)]
        elif len(indices) != len(cls.shape) - 1:
            raise IndexException("Rank mismatch.")
        return DataSet(cls, *indices)

    def _sympystr(self, p):
        """ For clarity the block number is printed"""
        return "%s_B%s" % (str(self.label), str(self.blocknumber))

    def simplelabel(self):
        """Returns the base label in case we need the label as a string to modify and create a new one with modifications
        """
        return "%s" % (str(self.label))

    @staticmethod
    def location():
        """ The location is the relative grid location, it is presently hard coded to the grid location.
        This is provided so that if in future staggered grid arrangement can be implemented
        """
        return [0 for i in range(DataSetBase.block.ndim)]


class DataSet(Indexed):
    """ It is the Data set of a DatasetBase which is defined on a block.

    :param base: name of the Dataset base, this is of type DataSetBase
    :param list indices: indices of the dataset, in out context it is the value of the dataset at the relative location from the dataset base location (example, [0,0] or [1,0])
    :returns: Dataset at the specified location
    :rtype: DataSet

    .. note::
        As given in SymPy's Indexed object documentation it is always advised to create a DataSetBase and then
        a Dataset
    """
    is_commutative = True
    is_Indexed = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *indices, **kwargs):
        if not isinstance(base, DataSetBase):
            raise ValueError("Declare DatasetBase and instantiate a dataset")
        ret = Indexed.__new__(cls, base, *indices)
        return ret

    def _sympystr(self, p):
        allinds = [i for i in self.indices if not isinstance(i, Idx)]
        indices = list(map(p.doprint, allinds))
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
        :rtype: list
        """
        return [i for i in self.indices if not isinstance(i, Idx)]


class GridIndex(Indexed):
    """ An Indexed object to get the local grid point
    """
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
    :rtype: CoordinateObject
    """

    def __new__(self, label, number, **assumptions):
        self._sanitize(assumptions, self)  # Remove any 'None's, etc.
        self.name = str(label) + str(number)
        # Make this into a new SymPy Symbol object.
        self = Symbol.__xnew__(self, self.name, **assumptions)
        self.number = number
        self.base = label
        return self

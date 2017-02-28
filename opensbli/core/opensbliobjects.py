
from sympy import Symbol, Eq, flatten, srepr
from sympy.tensor import Idx, IndexedBase, Indexed
from sympy import pprint
#class EinsteinIndex(Symbol):
    #return
class EinsteinTerm(Symbol):

    """ Represents any symbol in the equation as a SymPy Symbol object which in turn represents an Einstein term.
    This could be e.g. tau_i_j, but can also be e.g. u_i, rho.
    In other words, all symbols in the equation are Einstein terms, but they can have zero or more indices.
    Simple way of writing this class would be
    a. Notation, which gives the notaion as an indexed object u[i], tau[i,j] etc. etc.
    b. Value which is nothing but u_i will be u0, u1, u2
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
        self.rank = 0
        #pprint([self, type(self)])

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
        if len(self.get_indices()) >0:
            st = IndexedBase("%s"%str(self.get_base()))
            st = st[self.get_indices()]
            #st.is_commutative= False
            st.expression = self
            self.indexed_object = st
            return st
        else:
            st = self
            #st.is_commutative= False
            st.expression = st
            self.indexed_object = st
            return st
    def get_base(self):
        """ Return the base name.

        :returns: The base name.
        :rtype: str
        """
        return self.name.split('_')[0]

    #def __getitem__(self, idx, value):
        #return

    def apply_index(self, idx, value):
        if idx in self.get_indices():
            if self.get_base():
                newval = str(self).replace("_%s"%str(idx), str(value))
                val = type(self)(newval)
                # Store val index
                #if isinstance(self, CoordinateObject):
                val.direction = value
            else:
                val = int(value)
        else:
            val = self
        return val

    def apply_multiple_indices(self, indices, index_map):
        """
        This is used while having multiple indices and all are replaced
        """
        dictionary_indices = {}
        newval = str(self)
        for m in index_map:
            newval = newval.replace("_%s"%str(m[0]), str(m[1]))
            dictionary_indices[m[0]] = m[1]
        val = type(self)(newval)
        # Store val index
        if isinstance(self, CoordinateObject):
            if self.get_indices():
                val.direction = dictionary_indices[self.get_indices()[0]]
        return val
class Constant(object):
    pass

class ConstantObject(EinsteinTerm, Constant):
    is_commutative = True
    def __new__(cls, label, **kwargs):
        ret = super(ConstantObject,cls).__new__(cls, label, **kwargs)
        ret.is_constant = True
        return ret

class CoordinateObject(EinsteinTerm):
    is_commutative = True
    is_Atom = True
    def __new__(cls, label, **kwargs):
        ret = super(CoordinateObject,cls).__new__(cls, label)
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
    def __new__(cls, label, **kwargs):
        ret = super(MetricObject,cls).__new__(cls, label)
        if 'time' in kwargs:
            ret.timecoordinate = True
        else:
            ret.timecoordinate = False
        return ret

class DataObject(EinsteinTerm):
    """ This represents the objects this constitutes one of the basic terms if the OpenSBLI objects"""
    is_commutative = True
    #is_commutative = True
    is_Atom = True
    def __new__(cls, label, **kw_args):
        ret = super(DataObject, cls).__new__(cls, label, **kw_args)
        #ret.location = [0]*cls.ndim
        return ret
from sympy.core import Expr, Tuple, Symbol, sympify, S
from sympy.core.compatibility import is_sequence, string_types, NotIterable, range
from sympy.core.cache import cacheit
class DataSetBase(IndexedBase):
    is_commutative = True
    #is_Symbol = True
    #is_symbol = True
    is_Atom = True
    block = None
    @cacheit
    def __new__(cls, label, **kw_args):
        if not cls.block:
            raise ValueError("Set the block for DataSetBase")
        #sym = Symbol("%s_B%d"%(label, cls.block.blocknumber))
        sym = label
        shape = list(cls.block.shape) + [Idx(cls.block.blockname)]
        #cls.label = label
        ret = super(DataSetBase, cls).__new__(cls, sym, shape, **kw_args)
        ret.noblockname = Symbol(str(label))
        ret.blockname = cls.block.blockname
        ret.blocknumber = cls.block.blocknumber
        #ret.HDF5out = False
        #ret.HDF5inp = False
        return ret
    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h
    def _hashable_content(self):
        #print "IN ", self._args
        #exit()
        return str(self.label) + self.blockname
    def __getitem__(cls, indices, **kw_args):
        if len(indices) == len(cls.shape) and indices[-1] == Idx(cls.blockname):
            pass
        elif len(indices) == len(cls.shape) -1:
            indices += [Idx(cls.blockname)]
        elif len(indices) != len(cls.shape) -1:
            raise IndexException("Rank mismatch.")
        return DataSet(cls, *indices)
    def _sympystr(self, p):
        """ For clarity the block number is printed"""
        return "%s_B%s"%(str(self.label), str(self.blocknumber))
    def simplelabel(self):
        """Returns the abse label in case we need to add strings to the end of it
        e.g. Creating the old variables in Runge Kutta scheme
        """
        return "%s"%(str(self.label))
    #__xnew_cached_ = staticmethod(
        #cacheit(__new_stage2__))
    @staticmethod
    def location():
        return [0 for i in range(DataSetBase.block.ndim)]
    #def __eq__(

    #def _pretty(self, printer):
        """Pretty Printing method. """
        #from sympy.printing.pretty.stringpict import prettyForm
        #pforms = [printer._print("_B"), printer._print(str(self.blocknumber))]
        ##expr = prettyForm(*printer.join("", pforms))
        ##pform = prettyForm(printer._print(self.label))
        ##b = pform.baseline
        ##pform.baseline = pform.height() - 0.5
        #pform =prettyForm(*pform.right(*printer.join("", pforms)))
        ##pform.baseline = b
        #return pform


class DataSet(Indexed):
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
        allinds = [i for i in self.indices if not isinstance(i,Idx)]
        indices = list(map(p.doprint, allinds))
        #return "%s[%s]" % (p.doprint(self.base), ", ".join(indices))
        return "%s" % (p.doprint(self.base))
    def _pretty(self, printer):
        """Pretty Printing method. """
        from sympy.printing.pretty.stringpict import prettyForm
        allinds = [i for i in self.indices if not isinstance(i,Idx)]
        pforms = [printer._print(a) for a in allinds]
        expr = prettyForm(*printer.join(",", pforms).parens(left="[", right="]"))
        expr = prettyForm(*expr.left(printer._print(self.base)))
        return expr
    @property
    def get_grid_indices(self):
        return [i for i in self.indices if not isinstance(i, Idx)]

class ConstantIndexed(Indexed, Constant):
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
        ret = Indexed.__new__(cls, base , *indices)
        ret.is_constant = True
        return ret
    @property
    def location(cls):
        return list(cls.args[1:])

class GridIndex(Indexed):
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
    def free_symbols(self):
        return self
class GridIndexedBase(IndexedBase):
    is_commutative = True
    is_Symbol = True
    is_symbol = True
    is_Atom = True
    def __new__(cls, label, block, **kw_args):
        sym = label
        shape = list(block.shape)
        #cls.label = label
        ret = super(GridIndexedBase, cls).__new__(cls, sym, [block.ndim], **kw_args) # Shape would be of size ndim
        ret.blockname = block.blockname
        ret.blocknumber = block.blocknumber
        return ret
    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h
    def _hashable_content(self):
        return str(self.label) + self.blockname
    def __getitem__(cls, indices, **kw_args):
        if isinstance(indices, int):
            indices = [indices]
        if len(indices) != len(cls.shape):
            raise IndexException("Rank mismatch.")
        return GridIndex(cls, *indices)
    def _sympystr(self, p):
        """ For clarity the block number is printed"""
        return "%s_B%s"%(str(self.label), str(self.blocknumber))

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
        #pprint([self, type(self)])

        # Extract the indices, which are always preceded by an underscore.
        indices = self.name.split('_')[1:]
        self.indices = [Idx(x) for x in indices]
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
    def __new__(cls, label):
        ret = super(ConstantObject,cls).__new__(cls, label)
        ret.is_constant = True
        return ret

class CoordinateObject(EinsteinTerm):
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
    def __new__(cls, label, **kw_args):
        ret = super(DataObject, cls).__new__(cls, label, **kw_args)
        #ret.location = [0]*cls.ndim
        return ret
class DataSetBase(IndexedBase):
    #def modify_ranges(cls, direction_range, direction):
        #cls.ranges[direction] = direction_range
        #return
    def __init__(self, label, **kwargs):
        IndexedBase.__init__(label)
        return
    def set_ranges(self, range_is):
        self._ranges = range_is
        return
    def modify_ranges(self, direction_range, direction):
        self._ranges[direction] = direction_range
        return
    #def get_location_dataset(cls, location):
        #print(Eq(cls, type(cls)(str(cls), **{'location': location})))
        #ret = type(cls)(str(cls), **{'location': location})
        #return ret
class DataSet(Indexed):
    dimensions = 0
    def __new__(cls, label, **kwargs):
        if 'location' in kwargs:
            location = kwargs['location']
        else:
            location = [0]*cls.dimensions
        base = DataSetBase(label)
        ret = Indexed.__new__(cls, base , *location)
        return ret
    @property
    def location(cls):
        return list(cls.args[1:])

    def get_location_dataset(cls, location): ### Used for updating indices for weighting of points
        base = cls.base
        ret = type(cls)(str(base), **{'location': location})
        return ret
    @property
    def requires(cls):
        return list(cls.atoms(IndexedBase))[0]

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
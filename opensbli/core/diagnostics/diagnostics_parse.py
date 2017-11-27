
from sympy import pprint, Equality
from sympy.core.function import AppliedUndef, Function
#from opensbli.core.opensbliequations import OpenSBLIEquation as Eq
from opensbli.core.opensbliobjects import EinsteinTerm, DataObject
from opensbli.core.opensblifunctions import BasicDiscretisation

class ReductionVariable(EinsteinTerm, BasicDiscretisation):
    pass

class ReductionSum(ReductionVariable):
    def __new__(cls, name):
        ret = super(ReductionSum, cls).__new__(cls, name)
        return ret

class SpatialSum(Function, BasicDiscretisation):
    def __new__(cls, *args):
        ret = super(SpatialSum, cls).__new__(cls, *args)
        return ret

class SpatialMean(Function, BasicDiscretisation):
    def __new__(cls, *args):
        # check if the argument is of type spatial sum then donot apply again
        if isinstance(args[0], SpatialSum):
            a = args[0]
        else:
            a = SpatialSum(*args)
        ret = super(SpatialMean, cls).__new__(cls, a)
        return ret
    
    def get_lhs_variable(self, name):
        return ReductionSum(str(name))
    
    def discretise_spatial(self, block):
        
        return


class TemporalSum(Function, BasicDiscretisation):
    def __new__(cls, *args):
        ret = super(TemporalSum, cls).__new__(cls, *args)
        return ret

class TemporalMean(Function, BasicDiscretisation):
    def __new__(cls, *args):
        # check if the argument is of type spatial sum then donot apply again
        if isinstance(args[0], TemporalSum):
            a = args[0]
        else:
            a = TemporalSum(*args)
        ret = super(TemporalMean, cls).__new__(cls, a)
        return ret
    
    def discretise_spatial(self, block):
        
        return
    
    def discretise_temporal(self, block):
        
        return
    
    def get_lhs_variable(self, name):
        return DataObject(str(name))

class Mean(object):
    
    def __new__(cls, *args, **kwargs):
        if 'temporal' in kwargs:
            ret = TemporalMean(*args)
            #raise NotImplementedError("")
        else:
            ret = SpatialMean(*args)
        return ret
    
    def _sympystr(self, p):
        return "Mean(%s)" %( ','.join(str(a) for a in args))
    
class Diagnostic(Equality):
    
    def __new__(cls, *args, **kwargs):
        if (len(args) != 2):
            raise ValueError("")
        lhs = args[1].get_lhs_variable(args[0])
        rhs = args[1]
        ret = super(Diagnostic, cls).__new__(cls, lhs, rhs)
        return ret


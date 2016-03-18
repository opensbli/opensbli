from sympy import *
from sympy.parsing.sympy_parser import *

from .equations import *

class Kernel(object):
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
        self.inputoutput = {}
        #self.constants = {}
        self.classify_grid_objects()
        return
    def classify_grid_objects(self):
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
            self.inputoutput[v] = indexes
        idxs = flatten([list(eq.rhs.atoms(Idx)) for eq in self.equations])
        self.Idx = False
        if idxs:
            self.Idx = True
        self.constants = set(consts)
        return

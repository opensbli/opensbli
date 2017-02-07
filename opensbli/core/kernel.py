from sympy import flatten, Max
from .latex import *
from .opensbliobjects import DataSetBase, DataSet, ConstantIndexed, ConstantObject
def dataset_attributes(dset):
    """
    Move to datasetbase? Should we??
    """
    dset.block_number = None
    dset.read_from_hdf5 = False
    dset.dtype = None
    dset.size  = None
    dset.halo_ranges = None
    dset.block_name = None
    return dset

def constant_attributes(const):
    const.is_input = True
    const.dtype = None
    const.value = None
    return const

def copy_block_attributes(block, otherclass):
    """
    Move this to block
    """
    otherclass.block_number = block.blocknumber
    otherclass.ndim = block.ndim
    otherclass.block_name = block.blockname
    return
    ## similar define attributes for constant objects

class StencilObject(object):
    def __init__(self, name, stencil):
        self.name = name
        self.stencil = stencil
        return

class Kernel(object):

    """ A computational kernel which will be executed over all the grid points and in parallel. """
    mulfactor = {0:1, 1:1}
    opsc_access = {'ins':"OPS_READ", "outs": "OPS_WRITE", "inouts":"OPS_RW"}
    def __init__(self, block, computation_name = None):
        """ Set up the computational kernel"""
        copy_block_attributes(block,self)
        self.computation_name = computation_name
        self.kernel_no = block.kernel_counter
        self.kernelname = self.block_name + "Kernel%d"%self.kernel_no
        block.increase_kernel_counter
        self.equations = []
        self.halo_ranges = [[set(), set()] for d in range(block.ndim)]
        # self.stencil_names = {}
        return

    def set_computation_name(self, name):
        self.computation_name = name
        return

    def add_equation(self,equation):
        if isinstance(equation, list):
            self.equations += flatten([equation])
        elif isinstance(equation, Equality):
            self.equations += [equation]
        else:
            raise ValueError("Error in kernel add equation.")
        return

    def set_grid_range(self, block):
        self.ranges = block.ranges
        return

    #def get_max_halos(self, direction, side, block):
        #halos = block.boundary_halos[direction][side]
        #print halos
        #total_halos = 0
        #for h in halos:
            #total_halos = Max(total_halos, h.get_halos(side))
        #total_halos = self.mulfactor[side]*total_halos
        #return total_halos

    def get_max_halos(self,direction, side, block):
        halos =  block.boundary_halos
        halo_m = []
        halo_p = []
        for direction in range(len(halos)):
            max_halo_direction = []
            if halos[direction][0]:
                hal = [d.get_halos(0) for d in halos[direction][0]]
                halo_m += [min(hal)]
            else:
                halo_m += [0]
            if halos[direction][1]:
                hal = [d.get_halos(1) for d in halos[direction][1]]
                halo_p += [max(hal)]
            else:
                halo_p += [0]
        if side == 0:
            return halo_m[direction]
        elif side == 1:
            return halo_p[direction]

    def get_plane_halos(self, block):
        plane_halos = []
        for no, d in enumerate(block.boundary_halos):
            direction_halos = []
            direction_halos += [self.get_max_halos(no, 0, block)]
            direction_halos += [self.get_max_halos(no, 1, block)]
            plane_halos += [direction_halos]
        return plane_halos

    def set_boundary_plane_range(self, block, direction, side):
        self.ranges = block.ranges[:]
        self.ranges[direction] = block.ranges[direction][side]
        return

    def set_range(self, ranges):
        self.ranges = ranges # NOT REQUIRED
        return

    def set_halo_range(self, direction, side, types):
        #if not self.halo_ranges[direction][side]:
            #self.halo_ranges[direction][side] = set([types])
        #else:
        self.halo_ranges[direction][side].add(types)
        return

    def merge_halo_range(self, halo_range):
        for direction in range(len(self.halo_ranges)):
            self.halo_ranges[direction][0] = self.halo_ranges[direction][0] | halo_range[direction][0]
            self.halo_ranges[direction][1] = self.halo_ranges[direction][1] | halo_range[direction][1]

        return

    def set_kernel_evaluation_number(self, number, block):
        """
        This sets the evaluation number for the kernel so that they can be organised in
        algorithm
        """

        return

    def check_and_merge_kernels(self, kernel):
        """
        We donot check the equations only halo range is checked and updated
        """
        return
    @property
    def required_data_sets(self):
        requires = []
        for eq in self.equations:
            if isinstance(eq, Equality):
                requires += list(eq.rhs.atoms(DataSet))
        return requires
    @property
    def lhs_datasets(self):
        datasets = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                datasets = datasets.union(eq.lhs.atoms(DataSetBase))
        return datasets
    @property
    def rhs_datasets(self):
        datasets = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                datasets = datasets.union(eq.rhs.atoms(DataSetBase))
        return datasets
    @property
    def Rational_constants(self):
        rcs = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                rcs = rcs.union(eq.atoms(Rational))
        out = set()
        # Integers are also being returned as Rational numbers, remove any integers
        for rc in rcs:
            if not isinstance(rc, Integer):
                out.add(rc)
        return out
    @property
    def Inverse_constants(self):
        from sympy.core.function import _coeff_isneg
        # Only negative powers i.e. they correspond to division and they are stored into constant arrays
        inverse_terms = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                for at in eq.atoms(Pow):
                    if _coeff_isneg(at.exp) and not at.base.atoms(Indexed):
                        inverse_terms.add(at)
        return inverse_terms
    @property
    def constants(self):
        consts = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                consts = consts.union(eq.atoms(ConstantObject))
        return consts
    @property
    def IndexedConstants(self):
        consts = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                consts = consts.union(eq.atoms(ConstantIndexed))
        return consts

    def get_stencils(self):
        """ Returns the stencils for the datasets used in the kernel
        """
        stencil_dictionary = {}
        datasetbases = self.lhs_datasets.union(self.rhs_datasets)
        datasets = set()
        for eq in self.equations:
            if isinstance(eq, Equality):
                datasets = datasets.union(eq.atoms(DataSet))
        for s in datasets:
            if s.base in stencil_dictionary.keys():
                stencil_dictionary[s.base].add(tuple(s.indices))
            else:
                stencil_dictionary[s.base] = set()
                stencil_dictionary[s.base].add(tuple(s.indices))
        for key, val in stencil_dictionary.iteritems():
            stencil_dictionary[key] = frozenset(val)
        return stencil_dictionary

    def write_latex(self, latex):
        latex.write_string('The kernel is %s'%self.computation_name)
        #latex.write_string('. The range of evaluation is  %s \\ \n\n the halo ranges are %s'%(self.ranges, self.halo_ranges))
        for index, eq in enumerate(self.equations):
            if isinstance(eq, Equality):
                latex.write_expression(eq)
        return
    @property
    def opsc_code(self):
        block_name = self.block_name
        name = self.kernelname
        ins = self.rhs_datasets
        outs = self.lhs_datasets
        inouts = ins.intersection(outs)
        ins = ins.difference(inouts)
        outs = outs.difference(inouts)
        iter_range = "Testing"
        #print self.computation_name
        #pprint(self.stencil_names)
        code = ['ops_par_loop(%s, \"%s\", %s, %s, %s' % (name, self.computation_name, block_name, self.ndim, iter_range)]
        for i in ins:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)'%(i, 1, self.stencil_names[i], "double", self.opsc_access['ins'])]
        for o in outs:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)'%(o, 1, self.stencil_names[o], "double", self.opsc_access['outs'])]
        for io in inouts:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)'%(io, 1, self.stencil_names[io], "double", self.opsc_access['inouts'])]
        if self.IndexedConstants:
            for c in self.IndexedConstants:
                code += ["ops_arg_gbl(&%s, %d, \"%s\", %s)"%(c, 1, "double", self.opsc_access['ins'])]
        code = [',\n'.join(code) + ');\n\n']
        return code


    def ops_argument_call(self, array, stencil, precision, access_type):
        template = 'ops_arg_dat(%s, %d, %s, \"%s\", %s)'
        return template % (array, 1, stencil, self.dtype, access_type)

    def update_block_datasets(self, block):
        """
        Check the following
        a. existing.block_number is same as kernel
        b. set the range to block shape
        c. Update the halo ranges (similar to how we update the halo ranges of a kernel)

        Apply the datasetbase attributes to the dataset and update the parameters
        dataset_attributes(d)
        1. d.block_numner to kernel block number
        2. d.size = block shape
        3. d.halo_ranges to kernel halo ranges
        """
        self.stencil_names = {}
        dsets = self.lhs_datasets.union(self.rhs_datasets)
        for d in dsets:
            if str(d) in block.block_datasets.keys():
                dset = block.block_datasets[str(d)]
                for direction in range(len(dset.halo_ranges)):
                    dset.halo_ranges[direction][0] = dset.halo_ranges[direction][0] | self.halo_ranges[direction][0]
                    dset.halo_ranges[direction][1] = dset.halo_ranges[direction][1] | self.halo_ranges[direction][1]
                block.block_datasets[str(d)] = dset
                if block.blocknumber != dset.block_number:
                    raise ValueError("Block number error")
                if block.shape != dset.size:
                    raise ValueError("Shape error")
            else:
                # Update dataset attributes
                d = dataset_attributes(d)
                d.size = block.shape
                d.block_number = block.blocknumber
                d.halo_ranges = self.halo_ranges
                d.block_name = block.blockname
                # Add dataset to block datasets
                block.block_datasets[str(d)] = d
        # Update rational constant attributes
        rational_constants = self.Rational_constants.union(self.Inverse_constants)
        if rational_constants:
            for rc in rational_constants:
                if rc not in block.Rational_constants.keys():
                    next_rc = block.get_next_rational_constant
                    next_rc = constant_attributes(next_rc)
                    next_rc.is_input = False
                    next_rc.value = rc
                    block.Rational_constants[rc] = next_rc

        constants = self.constants
        if constants:
            for c in constants:
                if str(c) not in block.constants:
                    const_obj = constant_attributes(c)
                    # pprint(const_obj.__dict__)
                    block.constants[str(c)] = const_obj
            # pprint(block.constants)
        stens = self.get_stencils()
        for dset, stencil in stens.iteritems():
            if stencil not in block.block_stencils.keys():
                name = 'stencil_%d_%d' % (block.blocknumber, len(block.block_stencils.keys()))

                block.block_stencils[stencil] = StencilObject(name, stencil)
            if dset not in self.stencil_names:
                self.stencil_names[dset] = block.block_stencils[stencil].name
            else:
                self.stencil_names[dset].add(block.block_stencils[stencil].name)

        # pprint(block.block_stencils)
        #print "\n"
        #pprint(self.stencil_names)
        # for key, value in self.stencil_names.iteritems():
        #     print key, value, stens[key], block.block_stencils[stens[key]].name
        # pprint([self.stencil_names])
        # pprint(stens)
        return
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
    dset.halos = None
    return dset

def constant_attributes(const):
    const.is_input = True
    const.dtype = None
    const.value = None
    return const

    ## similar define attributes for constant objects
class Kernel(object):

    """ A computational kernel which will be executed over all the grid points and in parallel. """
    mulfactor = {0:-1, 1:1}
    opsc_access = {'ins':"OPS_READ", "outs": "OPS_WRITE", "inouts":"OPS_RW"}
    def __init__(self, block, computation_name = None):
        """ Set up the computational kernel"""
        self.block_number = block.blocknumber
        self.ndim = block.ndim
        self.computation_name = computation_name
        self.kernel_no = 0 #WARNING: update
        self.equations = []
        self.halo_ranges = [[set(), set()] for d in range(block.ndim)]
        return

    def set_computation_name(self, name):
        self.computation_name = name
        return

    def add_equation(self,equation):
        self.equations += flatten([equation])
        return

    def set_grid_range(self, block):
        self.ranges = block.ranges
        return

    def get_max_halos(self, direction, side, block):
        halos = block.boundary_halos[direction][side]
        total_halos = 0
        for h in halos:
            total_halos = Max(total_halos, h.get_halos(side))
        total_halos = self.mulfactor[side]*total_halos
        return total_halos

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
    @property
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
        block_name = "OpensbliBlock%d"%self.block_number
        ins = self.rhs_datasets
        outs = self.lhs_datasets
        inouts = ins.intersection(outs)
        ins = ins.difference(inouts)
        outs = outs.difference(inouts)
        stens = self.get_stencils
        # print self.computation_name, "\n"
        # print stens
        unique_stencils = set()
        for stencil in stens.values():
            unique_stencils.add(stencil)
        # pprint(unique_stencils)
        self.stencil_names = self.create_stencil_names(unique_stencils)
        #pprint(self.stencil_names)
        # pprint(self.stencil_names)
        name = "OpensbliKernel_block%d_kernel%d"%(self.block_number, self.kernel_no)
        iter_range = "Testing"
        code = ['ops_par_loop(%s, \"%s\", %s, %s, %s' % (name, self.computation_name, block_name, self.ndim, iter_range)]
        for i in ins:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)'%(i, 1, self.stencil_name(i, stens), "double", self.opsc_access['ins'])]
        for o in outs:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)'%(o, 1, self.stencil_name(o, stens), "double", self.opsc_access['outs'])]
        for io in inouts:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)'%(io, 1, self.stencil_name(io, stens), "double", self.opsc_access['inouts'])]
        self.declare_OPS_stencils()
        return code

    def sort_stencil_indices(self, index_set):
        """ Helper function for relative_stencil. Sorts the relative stencil. """
        dim = len(list(index_set)[0])
        sorted_index_set = sorted(index_set, key=lambda tup: tuple(tup[i] for i in range(dim)))
        return sorted_index_set

    def create_stencil_names(self, stencils):
        names = {}
        base_name = 'stencil_'
        for position, stencil in enumerate(stencils):
            names[stencil] = 'stencil_%d_%d_%d' % (self.block_number, self.kernel_no , position)
        return names
    def declare_OPS_stencils(self):
        dtype = 'int'
        # sten_format = 'ops_stencil %%s = ops_decl_stencil(%%d,%%d,%%s,\"%%s\")%s' % (self.end_of_statement)
        OPS_stencils = []
        for stencil, name in self.stencil_names.iteritems():
            sorted_stencil = self.sort_stencil_indices(stencil)
            OPS_stencils += flatten(sorted_stencil)
        #pprint(OPS_stencils)
        return
    # def declare_stencils(self):
    #     """ Declare all the stencils used in the code. We do not differentiate between the stencils for each block.

    #     :returns: The OPSC code declaring the stencil.
    #     :rtype: str
    #     """

    #     code = ['%s Declare all the stencils used ' % (self.line_comment)]
    #     dtype_int = 'int'
    #     sten_format = 'ops_stencil %%s = ops_decl_stencil(%%d,%%d,%%s,\"%%s\")%s' % (self.end_of_statement)
    #     for key, value in self.stencil_dictionary.iteritems():
    #         count = len(key.split(',')) / self.ndim
    #         # 'value' is the name in the stencil's format
    #         code += [self.array(dtype_int, value + "_temp", [key])]
    #         code += [sten_format % (value, self.ndim, count, value + "_temp", key)]
    #     return code

    def stencil_name(self, arr, stencils):
        return self.stencil_names[stencils[arr]]

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
                # Add dataset to block datasets
                block.block_datasets[str(d)] = d
        # Update rational constant attributes
        rational_constants = self.Rational_constants.union(self.Inverse_constants)
        # pprint(self.constants)
        if rational_constants:
            for rc in rational_constants:
                if rc in block.Rational_constants.keys():
                    #print "in existing"
                    #print block.Rational_constants[rc].__dict__
                    pass
                else:
                    next_rc = block.get_next_rational_constant
                    next_rc = constant_attributes(next_rc)
                    next_rc.is_input = False
                    next_rc.value = rc
                    block.Rational_constants[rc] = next_rc
        pprint(block.Rational_constants)

        # Update constant attributes
        # constants = self.constants
        # if constants:
        #     for c in constants:
        #         if str(c) in block.constants
        #         const_obj = ConstantObject(c)
        #         const_obj = constant_attributes(const_obj)
        #         pprint(const_obj.__dict__)
        #         exit()
        return
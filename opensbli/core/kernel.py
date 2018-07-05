"""@brief
   @authors Satya Pramod Jammy, David J Lusher
   @contributors
   @details
"""

from sympy import flatten, Equality, Indexed
from sympy import Rational, Pow, Integer
from opensbli.core.opensbliobjects import DataSet, ConstantIndexed, ConstantObject,\
    GlobalValue, GroupedPiecewise, Constant
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.grid import GridVariable, Grididx
from opensbli.core.datatypes import SimulationDataType
from sympy.core.function import _coeff_isneg
from opensbli.utilities.helperfunctions import get_min_max_halo_values, dataset_attributes
from opensbli.core.datatypes import Int
import copy

_known_equation_types = (GroupedPiecewise, OpenSBLIEq)


class ConstantsToDeclare(object):
    constants = []

    @staticmethod
    def add_constant(constant, value=None, dtype=None):
        """ Adds a constant or list of constants to be declared in the final program."""
        if isinstance(constant, Constant):
            if constant not in ConstantsToDeclare.constants:
                ConstantsToDeclare.constants += [constant]
        elif isinstance(constant, list):
            for c in constant:
                if c not in ConstantsToDeclare.constants:
                    ConstantsToDeclare.constants += [c]
        else:
            raise ValueError("Unknown type of constant")
        return


def copy_block_attributes(block, otherclass):
    """ Move this to block."""
    otherclass.block_number = block.blocknumber
    otherclass.ndim = block.ndim
    otherclass.block_name = block.blockname
    return


class StencilObject(object):
    def __init__(self, name, stencil, ndim):
        self.name = name
        self.stencil = stencil
        self.ndim = ndim
        self.dtype = Int()
        return

    def sort_stencil_indices(self):
        """ Helper function for relative_stencil, used in OPSC. Sorts the relative stencil locations."""
        index_set = self.stencil
        dim = len(list(index_set)[0])
        sorted_index_set = sorted(index_set, key=lambda tup: tuple(tup[i] for i in range(dim)))
        return sorted_index_set


class Kernel(object):
    """ A computational kernel, which will be executed over all the grid points in parallel."""
    mulfactor = {0: 1, 1: 1}
    opsc_access = {'ins': "OPS_READ", "outs": "OPS_WRITE", "inouts": "OPS_RW"}

    def __init__(self, block, computation_name=None):
        """ Set up the computational kernel"""
        copy_block_attributes(block, self)
        self.computation_name = computation_name
        self.kernel_no = block.kernel_counter
        self.kernelname = self.block_name + "Kernel%03d" % self.kernel_no
        block.increase_kernel_counter
        self.equations = []
        self.halo_ranges = [[set(), set()] for d in range(block.ndim)]
        return

    def set_computation_name(self, name):
        """ Sets the name of the computation for this kernel."""
        self.computation_name = name
        return

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h

    def _hashable_content(self):
        return str(self.kernelname)

    def add_equation(self, equation):
        """ Add an equation or list of equations to be evaluated inside this computational kernel."""
        if isinstance(equation, list):
            self.equations += flatten([equation])
        elif isinstance(equation, Equality):
            self.equations += [equation]
        elif isinstance(equation, GroupedPiecewise):
            self.equations += [equation]
        elif equation:
            pass
        else:
            raise ValueError("Error when adding equations to the kernel.")
        return

    def set_grid_range(self, block):
        """ Sets the kernel range equal to the block ranges."""
        self.ranges = copy.deepcopy(block.ranges)
        return

    def set_halo_range(self, direction, side, types):
        """ Sets the halo ranges for the kernel which extend beyond the grid range."""
        if isinstance(types, set):
            for s in types:
                self.halo_ranges[direction][side].add(s)
        else:
            self.halo_ranges[direction][side].add(types)
        return

    def merge_halo_range(self, halo_range):
        """ Merges the halo range for 2 kernels."""
        for direction in range(len(self.halo_ranges)):
            self.halo_ranges[direction][0] = self.halo_ranges[direction][0] | halo_range[direction][0]
            self.halo_ranges[direction][1] = self.halo_ranges[direction][1] | halo_range[direction][1]

        return

    @property
    def lhs_datasetbases(self):
        datasets = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                datasets = datasets.union(eq.lhs_datasetbases)
            elif isinstance(eq, Equality):
                print eq
                raise TypeError("Equality should be of types %s" % _known_equation_types)
        return datasets

    @property
    def rhs_datasetbases(self):
        datasets = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                datasets = datasets.union(eq.rhs_datasetbases)
            elif isinstance(eq, Equality):
                raise TypeError("Equality should be of types %s" % _known_equation_types)
        return datasets

    @property
    def Rational_constants(self):
        rcs = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                rcs = rcs.union(eq.atoms(Rational))
        out = set()
        # Integers are also being returned as Rational numbers, remove any integers
        for rc in rcs:
            if not isinstance(rc, Integer):
                out.add(rc)
        return out

    @property
    def Inverse_constants(self):
        # Only negative powers i.e. they correspond to division and they are stored into constant arrays
        inverse_terms = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                for at in eq.atoms(Pow):
                    if _coeff_isneg(at.exp) and not (at.base.atoms(Indexed) or isinstance(at, GridVariable)):
                        inverse_terms.add(at)
        return inverse_terms

    @property
    def constants(self):
        consts = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                consts = consts.union(eq.atoms(ConstantObject))
        return consts

    @property
    def IndexedConstants(self):
        consts = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                consts = consts.union(eq.atoms(ConstantIndexed))
        return consts

    @property
    def global_variables(self):
        globals_vars_lhs = set()
        globals_vars_rhs = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                globals_vars_lhs = globals_vars_lhs.union(eq.atoms(GlobalValue))
        return globals_vars_rhs, globals_vars_lhs

    @property
    def grid_indices_used(self):
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
                if eq.atoms(Grididx):
                    return True
        return False

    def get_stencils(self):
        """ Returns the stencils for the datasets used in the kernel."""
        stencil_dictionary = {}
        datasets = set()
        for eq in self.equations:
            if isinstance(eq, _known_equation_types):
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
        if isinstance(self.computation_name, str):
            name = self.computation_name
        else:
            name = latex.latexify_expression(self.computation_name, mode='inline')
        latex.write_string('The kernel is %s, block is %d' % (name, self.block_number))
        range_of_eval = self.total_range()
        latex.write_string('The ranges are %s' % (','.join([str(d) for d in flatten(range_of_eval)])))
        for index, eq in enumerate(self.equations):
            if isinstance(eq, Equality):
                latex.write_expression(eq)
            elif isinstance(eq, GroupedPiecewise):
                print "Should be doing latex for grouped piecewise"  # TODO
        return

    def total_range(self):
        range_of_eval = []
        if isinstance(self.halo_ranges, ConstantIndexed) and isinstance(self.ranges, ConstantIndexed):
            ranges = self.ranges.value_access_c
            halos = self.halo_ranges.value_access_c
            for a in list(zip(ranges, halos)):
                range_of_eval += [' + '.join(a)]
        elif isinstance(self.halo_ranges, ConstantIndexed) or isinstance(self.ranges, ConstantIndexed):
            raise NotImplementedError("handling ranges and halo_ranges of different types is not implemented")
        else:
            halo_m, halo_p = get_min_max_halo_values(self.halo_ranges)
            range_of_eval = [[0, 0] for r in range(self.ndim)]
            for d in range(self.ndim):
                range_of_eval[d][0] = self.ranges[d][0] + halo_m[d]
                range_of_eval[d][1] = self.ranges[d][1] + halo_p[d]
            range_of_eval = flatten(range_of_eval)
        return range_of_eval

    @property
    def opsc_code(self):
        """ Creates the OPSC code for a kernel."""
        block_name = self.block_name
        name = self.kernelname
        ins = self.rhs_datasetbases
        outs = self.lhs_datasetbases
        inouts = ins.intersection(outs)
        ins = ins.difference(inouts)
        outs = outs.difference(inouts)
        if len(self.equations) == 0:
            raise ValueError("Kernel %s does not have any equations." % self.computation_name)
        range_of_eval = self.total_range()
        dtype = Int().opsc()
        iter_name = "iteration_range_%d_block%d" % (self.kernel_no, self.block_number)
        iter_name_code = ['%s %s[] = {%s};' % (dtype, iter_name, ', '.join([str(s) for s in flatten(range_of_eval)]))]
        code = []
        # TODO check the dtype from the dataset
        sim_dtype = SimulationDataType.opsc()
        code += ['ops_par_loop(%s, \"%s\", %s, %s, %s' % (name, self.computation_name, block_name, self.ndim, iter_name)]
        for i in ins:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)' % (i, 1, self.stencil_names[i], sim_dtype, self.opsc_access['ins'])]  # WARNING dtype
        for o in outs:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)' % (o, 1, self.stencil_names[o], sim_dtype, self.opsc_access['outs'])]  # WARNING dtype
        for io in inouts:
            code += ['ops_arg_dat(%s, %d, %s, \"%s\", %s)' % (io, 1, self.stencil_names[io], sim_dtype, self.opsc_access['inouts'])]  # WARNING dtype
        if self.IndexedConstants:
            for c in self.IndexedConstants:
                code += ["ops_arg_gbl(&%s, %d, \"%s\", %s)" % (c, 1, sim_dtype, self.opsc_access['ins'])]
        if self.global_variables:
            # We need to write the size of an array for global indexed
            global_ins, global_outs = self.global_variables
            if global_ins.intersection(global_outs):
                raise NotImplementedError("Input output of global variables is not implemented")
            for c in global_ins:
                code += ["ops_arg_gbl(&%s, %d, \"%s\", %s)" % (c, 1, c.datatype.opsc(), self.opsc_access['ins'])]
            for c in global_outs:
                code += ["ops_arg_gbl(&%s, %d, \"%s\", %s)" % (c, 1, c.datatype.opsc(), self.opsc_access['outs'])]
        if self.grid_indices_used:
            code += ["ops_arg_idx()"]
        code = [',\n'.join(code) + ');\n\n']  # WARNING dtype
        code = iter_name_code + code
        return code

    def ops_argument_call(self, array, stencil, precision, access_type):
        template = 'ops_arg_dat(%s, %d, %s, \"%s\", %s)'
        return template % (array, 1, stencil, self.dtype, access_type)

    def update_block_datasets(self, block):
        """ Check the following
        a. existing.block_number is same as kernel
        b. set the range to block shape
        c. Update the halo ranges (similar to how we update the halo ranges of a kernel)

        Apply the datasetbase attributes to the dataset and update the parameters
        dataset_attributes(d)
        1. d.block_numner to kernel block number
        2. d.size = block shape
        3. d.halo_ranges to kernel halo ranges."""
        self.stencil_names = {}
        dsets = self.lhs_datasetbases.union(self.rhs_datasetbases)
        # New logic for the dataset delcarations across blocks
        for d in dsets:
            if str(d) in block.block_datasets.keys():
                dset = block.block_datasets[str(d)]
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
                d.halo_ranges = [[set(), set()] for d1 in range(block.ndim)]
                for direction in range(len(d.halo_ranges)):
                    d.halo_ranges[direction][0] = block.get_all_scheme_halos()
                    d.halo_ranges[direction][1] = block.get_all_scheme_halos()
                d.block_name = block.blockname
                # Add dataset to block datasets
                block.block_datasets[str(d)] = d

        stens = self.get_stencils()
        for dset, stencil in stens.iteritems():
            if stencil not in block.block_stencils.keys():
                name = 'stencil_%d_%02d' % (block.blocknumber, len(block.block_stencils.keys()))

                block.block_stencils[stencil] = StencilObject(name, stencil, block.ndim)
            if dset not in self.stencil_names:
                self.stencil_names[dset] = block.block_stencils[stencil].name
            else:
                self.stencil_names[dset].add(block.block_stencils[stencil].name)
        return

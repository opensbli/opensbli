"""@brief Algorithm generation
   @author Satya Pramod Jammy
   @contributors Satya Pramod Jammy and David J Lusher
   @details Impelments the tree-based structure with the attribute (components) controls the
   depth of a node.
"""

from sympy.core.compatibility import is_sequence
from sympy.printing.ccode import C99CodePrinter
from sympy.core.relational import Equality
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Constant, DataSetBase, DataObject, GroupedPiecewise
from opensbli.code_generation.algorithm.algorithm import Loop
from sympy import Symbol, flatten, pprint
from opensbli.core.grid import GridVariable
from opensbli.core.datatypes import SimulationDataType
from sympy import Pow, Idx

import os
import logging
LOG = logging.getLogger(__name__)
BUILD_DIR = os.getcwd()


class RationalCounter():

    # Counter for the kernels
    def __init__(self):
        self.name = 'rc%d'
        self.rational_counter = 0
        self.existing = {}

    @property
    def increase_rational_counter(self):
        self.rational_counter = self.rational_counter + 1
        return

    def get_next_rational_constant(self, numerical_value):
        from opensbli.core.kernel import ConstantsToDeclare
        name = self.name % self.rational_counter
        self.increase_rational_counter
        ret = ConstantObject(name)
        ret.value = numerical_value
        self.existing[numerical_value] = ret
        ConstantsToDeclare.add_constant(ret)
        return ret


rc = RationalCounter()


class OPSCCodePrinter(C99CodePrinter):

    """ Prints OPSC code. """
    dataset_accs_dictionary = {}
    settings_opsc = {'rational': False, 'kernel':False}

    def __init__(self, settings={}):
        """ Initialise the code printer. """
        self.settings_opsc = settings
        if 'rational' in settings.keys():
            self.settings_opsc = settings
        else:
            self.settings_opsc['rational'] = False
        C99CodePrinter.__init__(self, settings={})

    def _print_ReductionVariable(self, expr):
        return '*%s' % str(expr)

    def _print_Rational(self, expr):
        """Settings if rational is True then rational numbers are printed as they are.
        Else optimisations will be done for rational constants that are evaluated
        at the start of the program, to reduce divisions
        """
        if self.settings_opsc.get('rational', True):
            p, q = int(expr.p), int(expr.q)
            return '%d.0/%d.0' % (p, q)
        else:
            if expr in rc.existing:
                return self._print(rc.existing[expr])
            else:
                return self._print(rc.get_next_rational_constant(expr))

    def _print_Mod(self, expr):
        """All Mods are expressed as fmod currently and no iterger values"""
        args = map(ccode, expr.args)
        args = [x for x in args]
        result = ','.join(args)
        result = 'fmod(%s)' % result
        return result

    def _print_GridVariable(self, expr):
        """Prints the grid variable"""
        return str(expr)

    def _print_Max(self, expr):
        """MAXIMUM of the arguments, can handle any number of arguments:
        Max(a,b,c,d) is written as max(a, max(max(b,c),d))"""

        nargs = len(expr.args)
        args_code = [self._print(a) for a in expr.args]
        for i in range(nargs-1):
            # Max of the last 2 arguments in the array
            string_max = 'fmax(%s, %s)' % (args_code[-2], args_code[-1])
            # Remove the last 2 entries and append the max of the last 2
            del args_code[-2:]
            args_code.append(string_max)
        return str(args_code[0])

    def _print_DataObject(self, expr):
        """Raise error if a DataObject is found in the equation"""
        raise TypeError("Data object found in code generation, convert it to a dataset first, %s" % expr)
    
    def _print_Globalvariable(self, expr):
        # This should be handled in a different way if writing a kernel, if it is uses in the main cpp file it should be handled differently, only Global variables are supported and not arrays
        is_kernel = self.settings_opsc.get('kernel', False)
        if is_kernel:
            return ' *%s' %(str(expr))
        else:
            return "%s" % (str(expr))

    def _print_DataSetBase(self, expr):
        return str(expr)

    def _print_Equality(self, expr):
        return "%s == %s" % (self._print(expr.lhs), self._print(expr.rhs))

    def _print_DataSet(self, expr):
        """Prints the OpenSBLI dataset in the OPS format with the access numbers provided.
        Access numbers are updated for each kernel, see writing kernel in the OPSC class
        """
        base = expr.base
        if self.dataset_accs_dictionary[base]:
            indices = expr.get_grid_indices
            out = "%s[%s(%s)]" % (self._print(base), self.dataset_accs_dictionary[base].name, ','.join([self._print(i) for i in indices]))
            return out
        else:
            raise ValueError("Did not find the OPS Access for %s " % expr.base)

    def _print_Indexed(self, expr):
        """ Print out an Indexed object.

        :arg expr: The Indexed expression.
        :returns: The indexed expression, as OPSC code.
        :rtype: str
        """
        indices = [ind for ind in expr.indices]
        for number, index in enumerate(indices):
            for sym in index.atoms(Symbol):
                indices[number] = indices[number].subs({sym: 0})
        out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(index) for index in indices]))
        return out

    def _print_Grididx(self, expr):
        out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(expr.number)]))
        return out
    
    def _print_UserFunction(self, fn):
        """Prints the user defined funtion, we donot check if the funciton exists in the header files
        This allows for users to write their own OPS functions and use them in the code."""
        return self._print(fn.args[0]) + '(%s)' % (', '.join([self._print(arg) for arg in fn.args[1:]]))


def pow_to_constant(expr):
    """Finds the negative powers of constant objects and evaluates them to a new constant object
    to reduce divisions
    """

    from sympy.core.function import _coeff_isneg
    # Only negative powers i.e they correspond to division and they are stored into constant arrays
    inverse_terms = {}
    # change the name of Rational counter to rcinv
    orig_name = rc.name
    rc.name = 'rcinv%d'

    for at in expr.atoms(Pow):
        if _coeff_isneg(at.exp) and isinstance(at.base, ConstantObject):
            if at in rc.existing:
                inverse_terms[at] = rc.existing[at]
            else:
                inverse_terms[at] = rc.get_next_rational_constant(at)
    expr = expr.subs(inverse_terms)
    rc.name = orig_name  # change it back to the original name of Rational counter
    return expr


def ccode(expr, settings={}):
    """ Create an OPSC code printer object and write out the expression as an OPSC code string.

    :arg expr: The expression to translate into OPSC code.
    :arg Indexed_accs: Indexed OPS_ACC accesses.
    :arg constants: Constants that should be defined at the top of the OPSC code.
    :returns: The expression in OPSC code.
    :rtype: str
    """
    if isinstance(expr, Equality):
        if 'rational' in settings.keys():
            pass
        else:
            expr = pow_to_constant(expr)
        code_print = OPSCCodePrinter(settings)
        code = code_print.doprint(expr.lhs) \
            + ' = ' + OPSCCodePrinter(settings).doprint(expr.rhs)
        if isinstance(expr.lhs, GridVariable):
            code = code
        return code
    else:
        return OPSCCodePrinter(settings).doprint(expr)


class WriteString(object):
    def __init__(self, string):
        if isinstance(string, list):
            self.components = string
        elif isinstance(string, str):
            self.components = [string]
        else:
            raise ValueError("")
        return

    def __str__(self):
        return '\n'.join(self.components)

    def _print(self):
        return str(self)

    @property
    def opsc_code(self):
        return ['\n'.join(self.components)]


class OPSAccess(object):
    def __init__(self, no):
        self.name = "OPS_ACC%d" % no
        return


def indent_code(code_lines):
    """ Indent the code.

    :arg code_lines: The string or list of strings of lines of code to indent.
    :returns: A list of the indented line(s) of code.
    :rtype: list
    """

    p = C99CodePrinter()
    return p.indent_code(code_lines)


class OPSC(object):

    ops_headers = {'input': "const %s *%s", 'output': '%s *%s', 'inout': '%s *%s'}

    def __init__(self, algorithm):
        """ Generating an OPSC code from the algorithm"""
        if not algorithm.MultiBlock:
            self.MultiBlock = False
            self.dtype = algorithm.dtype
            # First write the kernels, with this we will have the Rational constants to declare
            self.write_kernels(algorithm)
            def_decs = self.opsc_def_decs(algorithm)
            end = self.ops_exit()
            algorithm.prg.components = def_decs + algorithm.prg.components + end
            code = algorithm.prg.opsc_code
            code = self.before_main(algorithm) + code
            self.name = 'taylor_green_vortex'
            f = open('opensbli.cpp', 'w')
            # code = indent_code(code)
            f.write('\n'.join(code))
            f.close()
            print "Successfully generated OPS C code"
        return

    def wrap_long_lines(self, code_lines):
        """ """
        limit = 120
        space = ' '
        formatted_code = []
        for code in code_lines:
            if len(code) > limit:
                # count the leading spaces so that we can use it at the end
                Leading_spaces = len(code) - len(code.lstrip())
                codelstrip = code.lstrip()
                length_ofstring = Leading_spaces
                split_lines = []
                # start the string with the leading spaces
                string = Leading_spaces*space
                for s in codelstrip.split(" "):
                    length_ofstring = len(string) + len(s)
                    if length_ofstring >= limit:
                        split_lines += [string]
                        # Increase the indentation by 4 spaces
                        string = (Leading_spaces + Leading_spaces)*space + s
                    else:
                        string = string + " " + s
                formatted_code += split_lines + [string]
                # Check the correctness of the code
                s1 = ''.join(split_lines + [string])
                s1 = "".join(s1.split())
                s2 = "".join(code.split())
                if s1 == s2:
                    pass
                else:
                    print(code)
                    print('\n'.join(split_lines + [string]))
                    raise ValueError("The code and the formatted line are not same")
            else:
                formatted_code += [code]

        return formatted_code

    def kernel_header(self, tuple_list):
        code = []
        dtype = SimulationDataType.opsc()
        for key, val in (tuple_list):
            # if any of the list has the datatype then use the data type
            if hasattr(key, "datatype") and key.datatype:
                code += [self.ops_headers[val] % (key.datatype.opsc(), key)]
            else:
                code += [self.ops_headers[val] % (dtype, key)]
        code = ', '.join(code)
        return code

    def kernel_computation_opsc(self, kernel):
        ins = kernel.rhs_datasetbases
        outs = kernel.lhs_datasetbases
        inouts = ins.intersection(outs)
        ins = ins.difference(inouts)
        outs = outs.difference(inouts)
        # eqs = kernel.equations
        all_dataset_inps = list(ins) + list(outs) + list(inouts)
        all_dataset_types = ['input' for i in ins] + ['output' for o in outs] + ['inout' for io in inouts]
        # add the global variables to the inputs and outputs
        global_ins, global_outs = kernel.global_variables
        if global_ins.intersection(global_outs):
            raise NotImplementedError("Input output of global variables is not implemented")
        all_dataset_inps += list(global_ins) + list(global_outs)
        all_dataset_types += ['input' for i in global_ins] + ['output' for o in global_outs] 
        # Use list of tuples as dictionary messes the order
        header_dictionary = zip(all_dataset_inps, all_dataset_types)
        if kernel.IndexedConstants:
            for i in kernel.IndexedConstants:
                header_dictionary += [tuple([(i.base), 'input'])]
        other_inputs = ""
        if kernel.grid_indices_used:
            # print kernel.grid_index_name
            other_inputs += ", const int *idx"  # WARNING hard coded here
        else:
            other_inputs = ''
        # print header_dictionary
        code = ["void %s(" % kernel.kernelname + self.kernel_header(header_dictionary) + other_inputs + ')' + '\n{']
        ops_accs = [OPSAccess(no) for no in range(len(all_dataset_inps))]
        OPSCCodePrinter.dataset_accs_dictionary = dict(zip(all_dataset_inps, ops_accs))
        # Find all the grid variables and declare them at the top 
        gridvariables = set()
        out = []
        for eq in kernel.equations:
            gridvariables = gridvariables.union(eq.atoms(GridVariable))
            if isinstance(eq, Equality):
                out += [ccode(eq, settings={'kernel': True}) + ';\n']
            elif isinstance(eq, GroupedPiecewise):
                for i, (expr, condition) in enumerate(eq.args):
                    if i == 0:
                        out += ['if (%s)' % ccode(condition, settings={'kernel': True}) + '{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings={'kernel': True}) + ';\n']
                        else:
                            out += [ccode(expr, settings={'kernel': True}) + ';\n']
                        out += ['}\n']
                    elif condition != True:
                        out += ['else if (%s)' % ccode(condition, settings={'kernel': True}) + '{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings={'kernel': True}) + ';\n']
                        else:
                            out += [ccode(expr, settings={'kernel': True}) + ';\n']
                        out += ['}\n']
                    else:
                        out += ['else{\n']
                        if is_sequence(expr):
                            for eqn in expr:
                                out += [ccode(eqn, settings={'kernel': True}) + ';\n']
                        else:
                            out += [ccode(expr, settings={'kernel': True}) + ';\n']
                        out += ['}\n']
            else:
                raise TypeError("Unclassified type of equation.")
        for gv in gridvariables:
            code += ["%s %s = 0.0;" %(SimulationDataType.opsc(), str(gv))]
        code += out + ['}'] # close Kernel
        OPSCCodePrinter.dataset_accs_dictionary = {}
        return code

    def write_kernels(self, algorithm):
        from opensbli.core.kernel import Kernel
        kernels = self.loop_alg(algorithm, Kernel)
        files = [open('%s_kernels.h' % b.block_name, 'w') for b in algorithm.block_descriptions]
        for f in files:
            name = ('%s_kernel_H' % b.block_name).upper()
            f.write('#ifndef %s\n' % name)
            f.write('#define %s\n' % name)
        for k in kernels:
            out = self.kernel_computation_opsc(k) + ['\n']
            out = indent_code(out)
            out = self.wrap_long_lines(out)
            files[k.block_number].write('\n'.join(out))
        for f in files:
            f.write("#endif\n")
        files = [f.close() for f in files]
        return

    def ops_exit(self):
        return [WriteString("ops_exit();")]

    def before_main(self, algorithm):
        out = ['#include <stdlib.h> \n#include <string.h> \n#include <math.h>']
        from opensbli.core.kernel import ConstantsToDeclare
        for d in ConstantsToDeclare.constants:
            if isinstance(d, ConstantObject):
                out += ["%s %s;" % (d.datatype.opsc(), d)]
            elif isinstance(d, ConstantIndexed):
                if not d.inline_array:
                    indices = ''
                    for s in d.shape:
                        indices = indices + '[%d]' % s
                    out += ["%s %s%s;" % (d.datatype.opsc(), d.base.label, indices)]
        for b in algorithm.block_descriptions:
            out += ['#define OPS_%dD' % b.ndim]
        out += ['#include \"ops_seq.h\"']
        for b in algorithm.block_descriptions:
            out += ['#include \"%s_kernels.h\"' % b.block_name]
        return out

    def opsc_def_decs(self, algorithm):
        defs = []
        decls = []
        # Add OPS_init to the declarations as it should be called before all ops
        decls += self.ops_init()
        # First process all the constants in the definitions
        from opensbli.core.kernel import ConstantsToDeclare
        for d in ConstantsToDeclare.constants:
            if isinstance(d, Constant):
                defs += self.define_constants(d)
                decls += self.declare_ops_constants(d)
        # Once the constants are done define and declare OPS dats
        output = defs + decls
        defs = []
        decls = []
        # Define and declare blocks
        for b in algorithm.block_descriptions:
            # output += self.define_block(b)
            output += self.declare_block(b)
        # output += defs + decls
        # Define and declare datasets on each block
        f = open('defdec_data_set.h', 'w')
        datasets_dec = []
        output += [WriteString("#include \"defdec_data_set.h\"")]
        for d in algorithm.defnitionsdeclarations.components:
            if isinstance(d, DataSetBase):
                datasets_dec += self.declare_dataset(d)
        f.write('\n'.join(flatten([d.opsc_code for d in datasets_dec])))
        f.close()
        output += [WriteString("// Define and declare stencils")]
        from opensbli.core.kernel import StencilObject
        for d in algorithm.defnitionsdeclarations.components:
            if isinstance(d, StencilObject):
                # pprint(d.__dict__)
                output += self.ops_stencils_declare(d)
        # Loop through algorithm components to include any halo exchanges
        from opensbli.core.bcs import Exchange
        exchange_list = self.loop_alg(algorithm, Exchange)
        if exchange_list:
            f = open('bc_exchanges.h', 'w')  # write BC_exchange code to a separate file
            exchange_code = []
            for e in exchange_list:
                call, code = self.bc_exchange_call_code(e)
                exchange_code += [code]
            f.write('\n'.join(flatten(exchange_code)))
            f.close()
            output += [WriteString("#include \"bc_exchanges.h\"")]  # Include statement in the code
        output += self.ops_partition()

        return output

    def ops_stencils_declare(self, s):
        out = []
        dtype = s.dtype.opsc()
        name = s.name + 'temp'
        sorted_stencil = s.sort_stencil_indices()
        out = [self.declare_inline_array(dtype, name, [st for st in flatten(sorted_stencil) if not isinstance(st, Idx)])]
        out += [WriteString('ops_stencil %s = ops_decl_stencil(%d,%d,%s,\"%s\");' % (s.name, s.ndim, len(s.stencil), name, name))]
        return out

    def ops_partition(self):
        """ Initialise an OPS partition for the purpose of multi-block and/or MPI partitioning.

        :returns: The partitioning code in OPSC format. Each line is a separate list element.
        :rtype: list
        """

        return [WriteString('// Init OPS partition'), WriteString('ops_partition(\"\");\n')]

    def ops_init(self, diagnostics_level=None):
        """ The default diagnostics level is 1, which offers no diagnostic information and should be used for production runs.
        Refer to OPS user manual for more information.

        :arg int diagnostics_level: The diagnostics level. If None, the diagnostic level defaults to 1.
        :returns: The call to ops_init.
        :rtype: list
        """
        out = [WriteString('// Initializing OPS ')]
        if diagnostics_level:
            self.ops_diagnostics = True
            return out + [WriteString('ops_init(argc,argv,%d);' % (diagnostics_level))]
        else:
            self.ops_diagnostics = False
            return out + [WriteString('ops_init(argc,argv,%d);' % (1))]

    def Exchange_code(self, e):
        # out =
        return

    def bc_exchange_call_code(self, instance):
        off = 0
        halo = 'halo'
        # instance.transfer_size = instance.transfer_from
        # Name of the halo exchange
        name = instance.name
        # self.halo_exchange_number = self.halo_exchange_number + 1
        code = ['// Boundary condition exchange code on %s direction %s %s' % (instance.block_name, instance.direction, instance.side)]
        code += ['ops_halo_group %s %s' % (name, ";")]
        code += ["{"]
        code += ['int halo_iter[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_size]), ";")]
        code += ['int from_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_from]), ";")]
        code += ['int to_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_to]), ";")]
        # dir in OPSC. WARNING: Not sure what it is, but 1 to ndim works.
        from_dir = [ind+1 for ind in range(len(instance.transfer_to))]
        to_dir = [ind+1 for ind in range(len(instance.transfer_to))]
        code += ['int from_dir[] = {%s}%s' % (', '.join([str(ind) for ind in from_dir]), ";")]
        code += ['int to_dir[] = {%s}%s' % (', '.join([str(ind) for ind in to_dir]), ";")]
        # Process the arrays
        for no, arr in enumerate(instance.transfer_arrays):
            from_array = instance.from_arrays[no]
            to_array = instance.to_arrays[no]
            code += ['ops_halo %s%d = ops_decl_halo(%s, %s, halo_iter, from_base, to_base, from_dir, to_dir)%s'
                     % (halo, off, from_array.base, to_array.base, ";")]
            off = off+1
        code += ['ops_halo grp[] = {%s}%s' % (','.join([str('%s%s' % (halo, of)) for of in range(off)]), ";")]
        code += ['%s = ops_decl_halo_group(%d,grp)%s' % (name, off, ";")]
        code += ["}"]
        # Finished OPS halo exchange, now get the call
        instance.call_name = 'ops_halo_transfer(%s)%s' % (name, ";")
        call = ['// Boundary condition exchange calls', 'ops_halo_transfer(%s)%s' % (name, ";")]
        for no, c in enumerate(code):
            code[no] = WriteString(c).opsc_code
        return call, code

    def loop_alg(self, algorithm, type_of_component):
        type_list = []

        def _generate(components, type_list):
            for component1 in components:
                if hasattr(component1, 'components'):
                    _generate(component1.components, type_list)
                elif isinstance(component1, type_of_component):
                    if component1 in type_list:
                        pass
                    else:
                        type_list += [component1]
        for c in algorithm.prg.components:
            _generate([c], type_list)
        return type_list

    def define_block(self, b):
        if not self.MultiBlock:
            return [WriteString("ops_block %s;" % b.block_name)]
        else:
            raise NotImplementedError("")

    def declare_block(self, b):
        if not self.MultiBlock:
            out = [WriteString("// Define and Declare OPS Block")]
            out += [WriteString('ops_block %s = ops_decl_block(%d, \"%s\");' % (b.block_name, b.ndim, b.block_name))]
            return out
        else:
            raise NotImplementedError("")

    def define_constants(self, c):
        # Fix spacing on constant declarations %s=%s
        if isinstance(c, ConstantObject):
            if not isinstance(c.value, str):
                return [WriteString("%s = %s;" % (str(c), ccode(c.value, settings={'rational': True})))]
            else:
                return [WriteString("%s=%s;" % (str(c), c.value))]
        elif isinstance(c, ConstantIndexed):
            out = []
            if c.value:
                if len(c.shape) == 1:
                    if c.inline_array:
                        values = [ccode(c.value[i], settings={'rational': True}) for i in range(c.shape[0])]
                        return [WriteString("%s %s[] = {%s};" % (c.datatype.opsc(), c.base.label, ', '.join(values)))]
                    else:
                        indices = ''
                        for s in c.shape:
                            indices = indices + '[%d]' % s
                        out += [WriteString("%s %s%s;" % (c.base.label, indices))]
                        for i in range(c.shape[0]):
                            out += [WriteString("%s[%d] = %s;" % (str(c.base.label), i, ccode(c.value[i], settings={'rational': True})))]
                        return out
                else:
                    raise NotImplementedError("Indexed constant declaration is done for only one ")
            else:
                raise NotImplementedError("")
        else:
            print c
            raise ValueError("")

    def declare_ops_constants(self, c):
        if isinstance(c, ConstantObject):
            return [WriteString("ops_decl_const(\"%s\" , 1, \"%s\", &%s);" % (str(c), c.datatype.opsc(), str(c)))]
        elif isinstance(c, ConstantIndexed):
            return []
        return

    def declare_inline_array(self, dtype, name, values):
        return WriteString('%s %s[] = {%s};' % (dtype, name, ', '.join([str(s) for s in values])))

    def update_inline_array(self, name, values):
        out = []
        for no, v in enumerate(values):
            out += [WriteString("%s[%d] = %s;" % (name, no, v))]
        return out

    def define_dataset(self, dset):
        if not self.MultiBlock:
            return [WriteString("ops_dat %s;" % (dset))]

    def get_max_halos(self, halos):
        halo_m = []
        halo_p = []
        for direction in range(len(halos)):
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
        return halo_m, halo_p

    def declare_dataset(self, dset):
        declaration = WriteString("ops_dat %s;" % dset)
        out = [declaration, WriteString("{")]

        if dset.dtype:
            dtype = dset.dtype
        else:
            dtype = SimulationDataType.dtype()

        if dset.read_from_hdf5:
            temp = '%s = ops_decl_dat_hdf5(%s, 1, \"%s\", \"%s\", \"%s\");' % (dset,
                                                                               dset.block_name, dtype.opsc(), dset, dset.input_file_name)
            out += [WriteString(temp)]
        else:
            hm, hp = self.get_max_halos(dset.halo_ranges)
            halo_p = self.declare_inline_array("int", "halo_p", hp)
            halo_m = self.declare_inline_array("int", "halo_m", hm)
            sizes = self.declare_inline_array("int", "size", [str(s) for s in (dset.size)])
            base = self.declare_inline_array("int", "base", [0 for i in range(len(dset.size))])
            value = WriteString("%s* value = NULL;" % dtype.opsc())
            temp = '%s = ops_decl_dat(%s, 1, size, base, halo_m, halo_p, value, \"%s\", \"%s\");' % (dset,
                                                                                                     dset.block_name, dtype.opsc(), dset)
            out += [halo_p, halo_m, sizes, base, value, WriteString(temp)]

        out += [WriteString("}")]
        return out

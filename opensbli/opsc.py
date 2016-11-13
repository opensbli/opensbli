#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

from sympy import *
from sympy.printing.ccode import CCodePrinter
import os
from string import Template
from .equations import EinsteinTerm
from .diagnostics import ReductionVariable
import logging
LOG = logging.getLogger(__name__)
BUILD_DIR = os.getcwd()

import subprocess

try:
    from ops_translator.c import ops as translator
    have_ops = True
    print "Found translator module: ", translator
except ImportError:
    logging.warning("Could not import the OPS library. The generated OPSC code will need to be manually put through the translator.")
    have_ops = False


class OPSCCodePrinter(CCodePrinter):

    """ Prints OPSC code. """

    def __init__(self, Indexed_accs, constants):
        """ Initialise the code printer. """

        settings = {}
        CCodePrinter.__init__(self, settings)

        # Indexed access numbers are required in dictionary
        self.Indexed_accs = Indexed_accs
        self.constants = constants

    def _print_ReductionVariable(self, expr):
        return '*%s' % str(expr)

    def _print_Rational(self, expr):
        if self.constants is not None:
            if expr in self.constants.keys():
                return self.constants[expr]
            else:
                con = 'rc%d' % len(self.constants.keys())
                self.constants[expr] = con
                return self.constants[expr]
        else:
            p, q = int(expr.p), int(expr.q)
            return '%d.0/%d.0' % (p, q)

    def _print_Mod(self, expr):
        args = map(ccode, expr.args)
        args = [x for x in args]
        result = ','.join(args)
        result = 'fmod(%s)' % result
        return result

    def _print_Indexed(self, expr):
        """ Print out an Indexed object.

        :arg expr: The Indexed expression.
        :returns: The indexed expression, as OPSC code.
        :rtype: str
        """

        # Replace the symbols in the indices that are not time with `zero'
        indices = [ind for ind in expr.indices if ind != EinsteinTerm('t')]
        for number, index in enumerate(indices):
            for sym in index.atoms(Symbol):
                indices[number] = indices[number].subs({sym: 0})
        if self.Indexed_accs:
            if self.Indexed_accs[expr.base]:
                out = "%s[%s(%s)]" % (self._print(expr.base.label), self.Indexed_accs[expr.base], ','.join([self._print(index) for index in indices]))
            else:
                out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(index) for index in indices]))
        else:
            out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(index) for index in indices]))
        return out


def pow_to_constant(expr, constants):
    from sympy.core.function import _coeff_isneg
    # Only negative powers i.e. they correspond to division and they are stored into constant arrays
    inverse_terms = {}
    for at in expr.atoms(Pow):
        if _coeff_isneg(at.exp) and not at.base.atoms(Indexed):
            if at not in constants.keys():
                constants[at] = 'rinv%d' % len(constants.keys())
            inverse_terms[at] = constants[at]
    expr = expr.subs(inverse_terms)
    return expr, constants


def ccode(expr, Indexed_accs=None, constants=None):
    """ Create an OPSC code printer object and write out the expression as an OPSC code string.

    :arg expr: The expression to translate into OPSC code.
    :arg Indexed_accs: Indexed OPS_ACC accesses.
    :arg constants: Constants that should be defined at the top of the OPSC code.
    :returns: The expression in OPSC code.
    :rtype: str
    """
    if isinstance(expr, Eq):
        if constants:
            expr, constants = pow_to_constant(expr, constants)
        code_print = OPSCCodePrinter(Indexed_accs, constants)
        code = code_print.doprint(expr.lhs) \
            + ' = ' + OPSCCodePrinter(Indexed_accs, constants).doprint(expr.rhs)
        return code, code_print.constants
    return OPSCCodePrinter(Indexed_accs, constants).doprint(expr)


class OPSC(object):

    """ A class describing the OPSC language, and various templates for OPSC code structures (e.g. loops, declarations, etc). """

    # OPS Access types, used for kernel call
    ops_access = {'inputs': 'OPS_READ', 'outputs': 'OPS_WRITE', 'inputoutput': 'OPS_RW', 'reduction': 'OPS_INC'}
    # OPS kernel headers
    ops_header = {'inputs': 'const %s *%s', 'outputs': '%s *%s', 'inputoutput': '%s *%s', 'Idx': 'const int *%s',
                  'reduction': '%s *%s'}
    # Single line comment
    line_comment = "//"
    # Block/multi-line comment
    block_comment = ['/*', '*/']
    # End of statement delimiter
    end_of_statement = ";"
    # Commonly used brackets
    left_brace = "{"
    right_brace = "}"
    left_parenthesis = "("
    right_parenthesis = ")"

    def __init__(self, grid, spatial_discretisation, temporal_discretisation, boundary_condition, initial_conditions, IO, simulation_parameters, diagnostics=None):
        self.check_consistency(grid, spatial_discretisation, temporal_discretisation, boundary_condition, initial_conditions, IO)
        if diagnostics:
            self.diagnostics = self.to_list(diagnostics)
        else:
            self.diagnostics = None
        self.simulation_parameters = simulation_parameters
        # Update the simulation parameters from that of the grid
        for g in self.grid:
            self.simulation_parameters.update(g.grid_data_dictionary)
        self.initialise_ops_parameters()
        self.template()
        if have_ops:
            self.translate()
        return

    def initialise_ops_parameters(self):
        """ This initialises various OPS parameters like the name of the computational files,
        computation kernel name, iteration range, stencil name, etc. Most of these are specific to OPS.
        """

        # Multiblock or single block?
        if len(self.grid) == 1:
            self.multiblock = False
            self.nblocks = 1
        else:
            self.multiblock = True
            self.nblocks = len(self.grid)

        # Dimensions of the blocks
        ndim = list(set([len(self.grid[i].shape) for i in range(self.nblocks)]))
        if len(ndim) != 1:
            raise ValueError("Mismatch in the grid shape of the blocks.")
        self.ndim = ndim[0]
        name = self.simulation_parameters["name"]

        # Name of the block for OPSC
        self.block_name = '%s_block' % name

        # File names for computational kernels, each block will have its own computational kernel file
        self.computational_routines_filename = ['%s_block_%d_kernel.h' % (name, block) for block in range(self.nblocks)]

        # Kernel names for each block. This will be the file name + the kernel number
        self.computational_kernel_names = ['%s_block%d_%%d_kernel' % (name, block) for block in range(self.nblocks)]

        # Name for exchange boundary condition, stencil, iteration range and kernel name
        self.iteration_range_name = 'iter_range%d'

        # Name for the commonly-used stencils
        self.stencil_name = 'stencil%d'
        self.stencil_number = 0
        self.iteration_range_index = 0
        self.kernel_name_number = [0 for block in range(self.nblocks)]

        # Name for exchange boundary conditions
        self.halo_exchange_number = 0
        self.halo_exchange_name = 'halo_exchange%d'

        # Grid based arrays used for declaration and definition in OPSC format
        self.grid_based_arrays = set()

        # The global constants that are to be declared.
        self.constants = set()

        # OPS constants. These are the constants in the above list to be defined in OPS format.
        self.constant_values = {}
        self.rational_constants = {}

        # Reduction variables
        self.reduction_variables = set()

        # Dictionary of stencils. The key will be a stencil, and the value is the name of stencil.
        self.stencil_dictionary = {}

        # Data type of arrays
        self.dtype = self.simulation_parameters['precision']

        # Create the code directory
        if not os.path.exists(BUILD_DIR+'/%s_opsc_code' % name):
            os.makedirs(BUILD_DIR+'/%s_opsc_code' % name)
        self.CODE_DIR = BUILD_DIR + '/%s_opsc_code' % name
        return

    def template(self):
        """ Define the algorithm in pseudo-code and get all the code. """

        OPS_template = """\
        $header
        \n$main_start
            \n$initialise_constants
            \n$ops_init
            \n$declare_ops_constants
            \n$define_block
            \n$initialise_block
            \n$define_dat
            \n$initialise_dat
            \n$declare_stencils
            \n$declare_reductions
            \n$bc_exchange
            \n$ops_partition
            \n$initialisation
            \n$bc_calls
            \n$timer_start
                \n$timeloop
                    \n$time_start_calls
                        \n$innerloop
                            \n$time_calls
                            \n$bc_calls
                        \n$end_inner_loop
                    \n$time_end_calls
                    \n$io_time
                \n$end_time_loop
            \n$timer_end
            \n$print_timings
            \n$io_calls
            \n$ops_exit
        \n$main_end
        """

        # OPS template is written in multiple lines for clarity. Now, remove all white spaces for nicety of file printing.
        OPS_template = OPS_template.replace(" ", "")

        # Dictionary to store evaluated things in the code_template
        code_dictionary = {}

        # Convert OPS_template to a Python Template
        code_template = Template(OPS_template)

        # Start populating the code dictionary with the corresponding code

        # Get the main start and main end code
        code_dictionary['main_start'] = '\n'.join(self.main_start())
        code_dictionary['main_end'] = self.right_brace

        # Get the ops_init, ops_exit (footer) calls
        code_dictionary['ops_init'] = '\n'.join(self.ops_init())
        code_dictionary['ops_exit'] = '\n'.join(self.footer())
        code_dictionary['ops_partition'] = '\n'.join(self.ops_partition())

        # Set up the timers
        timer = self.ops_timers()
        code_dictionary['timer_start'] = '\n'.join(timer[0])
        code_dictionary['timer_end'] = '\n'.join(timer[1])
        code_dictionary['print_timings'] = '\n'.join(timer[2])

        # Define the main time loop
        name = 'iteration'  # Name for the iteration.
        code_dictionary['timeloop'] = self.loop_open(name, (0, self.simulation_parameters['niter'])) + '\n'
        code_dictionary['end_time_loop'] = self.loop_close()

        # Declare and initialise OPS block
        code_dictionary['define_block'] = '\n'.join(self.define_block())
        code_dictionary['initialise_block'] = '\n'.join(self.initialise_block())

        # Set up the time-stepping sub loop, if required.
        if self.temporal_discretisation[0].nstages > 1:
            number_of_stages = self.temporal_discretisation[0].nstages
            stage = self.temporal_discretisation[0].scheme.stage
            code_dictionary['innerloop'] = self.loop_open(stage, (0, number_of_stages)) + '\n'
            code_dictionary['end_inner_loop'] = self.loop_close()
            # Update constant values
            coeffs = self.temporal_discretisation[0].scheme.get_coefficients()
            for key, value in coeffs.iteritems():
                self.simulation_parameters[str(key)] = value
        else:
            code_dictionary['innerloop'] = ""
            code_dictionary['end_inner_loop'] = ""

        # Get the computational routines
        computational_routines = self.get_block_computations()
        # Write the computational routines to block computation files
        self.write_computational_routines(computational_routines)

        # Computation calls
        # First the inner computation calls
        computations = [self.spatial_discretisation[block].computations + self.temporal_discretisation[block].computations for block in range(self.nblocks)]
        calls = self.get_block_computation_kernels(computations)
        code_dictionary['time_calls'] = '\n'.join(['\n'.join(calls[block]) for block in range(self.nblocks)])

        # Computations at the start of the time stepping loop
        computations = [self.temporal_discretisation[block].start_computations if self.temporal_discretisation[block].start_computations else [] for block in range(self.nblocks)]
        calls = self.get_block_computation_kernels(computations)
        code_dictionary['time_start_calls'] = '\n'.join(['\n'.join(calls[block]) for block in range(self.nblocks)])

        # Computations at the end of the time stepping loop
        computations = [self.temporal_discretisation[block].end_computations if self.temporal_discretisation[block].end_computations else [] for block in range(self.nblocks)]
        calls = self.get_block_computation_kernels(computations)
        code_dictionary['time_end_calls'] = '\n'.join(['\n'.join(calls[block]) for block in range(self.nblocks)])

        # computations for the initialisation, if there are computations, later this should be included to read from file
        computations = [self.initial_conditions[block].computations if self.initial_conditions[block].computations else [] for block in range(self.nblocks)]
        calls = self.get_block_computation_kernels(computations)
        code_dictionary['initialisation'] = '\n'.join(['\n'.join(calls[block]) for block in range(self.nblocks)])

        # Do the exchange boundary conditions
        code_dictionary = self.update_boundary_conditions(code_dictionary)

        # IO calls
        code_dictionary = self.get_io(code_dictionary)

        # Get diagnostics kernel calls
        if self.diagnostics:
            code_dictionary = self.get_diagnostic_kernels(code_dictionary)

        # Define and initialise all the data arrays used in the computations
        code_dictionary['define_dat'] = '\n'.join(self.define_dat())
        code_dictionary['initialise_dat'] = '\n'.join(self.initialise_dat())

        # Header
        code_dictionary['header'] = '\n'.join(self.header())

        # Stencils
        code_dictionary['declare_stencils'] = '\n'.join(self.declare_stencils())

        # Initialise constants and Declare  constants in OPS format
        code_dictionary['initialise_constants'] = '\n'.join(self.initialise_constants())
        code_dictionary['declare_ops_constants'] = '\n'.join(self.declare_ops_constants())

        # Reduction declarations
        code_dictionary['declare_reductions'] = '\n'.join(self.declare_reduction_variables())

        # Write the main file
        code_template = code_template.safe_substitute(code_dictionary)
        self.write_main_file(code_template)
        return

    def get_diagnostic_kernels(self, code_dictionary):
        """ Loop over blocks, loop over each diagnostics object (can be reduction etc.), and get the kernel call.
        If it is a Reduction, get the reduction result and write the output to a file. """

        from .diagnostics import Reduction as R
        for block in range(self.nblocks):
            for diagnostic in self.diagnostics[block]:
                if isinstance(diagnostic, R):
                    calls = []
                    for computation in diagnostic.computations:
                        calls += self.kernel_call(computation)
                        if computation.reductions:
                            calls += self.get_reduction_results(computation.reductions)
                            calls += self.print_reduction_results(computation.reductions)
                    if diagnostic.compute_every:
                        l = ccode(Mod('iteration', diagnostic.compute_every))
                        calls = ['if(%s == 0)' % l] + [self.left_brace] + calls
                        calls += [self.right_brace]
                        code_dictionary['io_time'] += '\n'.join(calls)
                    else:
                        code_dictionary['io_calls'] += '\n'.join(calls)

        return code_dictionary

    def print_reduction_results(self, reductions):
        """ Prints the reduction results, as these are called at the end of the simulation
        prints time + all reduction results in a single line"""

        template = "ops_printf(\"%s\\n\", %s)%s"
        all_reductions = '%g, ' + ', '.join([str('%g') for red in reductions])
        all_reduction_results = '(iteration + 1)*deltat, ' + ', '.join([str('%s_reduction') % red for red in reductions])
        return [template % (all_reductions, all_reduction_results, self.end_of_statement)]

    def get_reduction_results(self, reductions):
        """ Returns the code for OPS reduction result"""

        template = "%s %s_reduction = 0.0%s \n ops_reduction_result(%s, &%s_reduction)%s"
        return [template % (self.dtype, red, self.end_of_statement, red, red, self.end_of_statement)
                for red in reductions]

    def declare_reduction_variables(self):
        template = "ops_reduction %s = ops_decl_reduction_handle(sizeof(%s), \"%s\", \"reduction_%s\")%s"
        return [template % (red, self.dtype, self.dtype, red, self.end_of_statement) for red in self.reduction_variables]

    def initialise_constants(self):
        """ Initialise all constant values. """

        constant_initialisation = []
        constant_dictionary = {}
        values = [self.simulation_parameters[str(constant)] for constant in self.constants]
        constant_dictionary = dict(zip(self.constants, values))
        # Sort constants
        sorted_constants = []
        sorted_constants = self.sort_constants(constant_dictionary, sorted_constants)

        for constant in sorted_constants:
            val = self.simulation_parameters[str(constant)]
            if isinstance(constant, IndexedBase):
                if constant.ranges != len(val):
                    raise ValueError("The indexed constant %s should have only %d values" % (constant, constant.ranges))
                for r in range(constant.ranges):
                    constant_initialisation += ["%s[%d] = %s%s" % (constant, r, ccode(val[r]), self.end_of_statement)]
            else:
                constant_initialisation += ["%s = %s%s" % (constant, ccode(val), self.end_of_statement)]
        return constant_initialisation

    def sort_constants(self, constant_dictionary, sorted_constants):
        """ Sorts the constants. The function breaks , if unable to sort in 1000 iterations.
        Prints out various stuff"""
        types_known = (float, int, Rational)

        # Sort the Float, rational and integer constants
        sorted_constants = [key for key, value in constant_dictionary.iteritems()
                            if isinstance(value, types_known) and key not in sorted_constants]

        # Known types constants for the Indexed Constants
        sorted_constants += [key for key, value in constant_dictionary.iteritems()
                             if isinstance(key, IndexedBase) and all(isinstance(v, types_known) for v in value)
                             and key not in sorted_constants]

        # get the other constants that are symbol/ Einstein Term
        key_list = [key for key in constant_dictionary.keys() if key not in sorted_constants]

        requires_list = [constant_dictionary[key].atoms(Symbol) for key in key_list]
        zipped = zip(key_list, requires_list)

        iter_count = 0
        while key_list:
            iter_count = iter_count+1
            sorted_constants += [x for (x, y) in zipped if all(req in sorted_constants for req in y)]
            key_list = [key for key in constant_dictionary.keys() if key not in sorted_constants]
            requires_list = [constant_dictionary[key].atoms(Symbol) for key in key_list]
            zipped = zip(key_list, requires_list)
            if iter_count > 1000:
                pprint("Constant sorting recursion reached")
                pprint([req for req in requires_list[0]])
                pprint(key_list)
                pprint([(req in sorted_constants, req) for req in requires_list[0]])
                raise ValueError("Exiting sort evaluations ")
        return sorted_constants

    def declare_ops_constants(self):
        """ Declare each constant as an OPS constant.

        :returns: A list of OPS constant declarations.
        :rtype: list
        """

        ops_const = []
        for constant in self.constants:
            if not isinstance(constant, IndexedBase) and isinstance(constant, str):
                ops_const += ["ops_decl_const(\"%s\" , 1, \"%s\", &%s)%s" % (constant, self.dtype, constant, self.end_of_statement)]
            elif constant.is_integer:
                ops_const += ["ops_decl_const(\"%s\" , 1, \"int\", &%s)%s" % (constant, constant, self.end_of_statement)]
            elif not isinstance(constant, IndexedBase):
                ops_const += ["ops_decl_const(\"%s\" , 1, \"%s\", &%s)%s" % (constant, self.dtype, constant, self.end_of_statement)]
        return ops_const

    def write_main_file(self, code_template):
        """ Write the main .cpp file. The base name of the file will be the same as the simulation's name.

        :arg code_template: The templated code to write out.
        :returns: None
        """

        mainfile = open(self.CODE_DIR+'/'+'%s.cpp' % self.simulation_parameters["name"], 'w')
        code_template = self.indent_code(code_template)
        mainfile.write(code_template)
        mainfile.close()
        return

    def indent_code(self, code_lines):
        """ Indent the code.

        :arg code_lines: The string or list of strings of lines of code to indent.
        :returns: A list of the indented line(s) of code.
        :rtype: list
        """

        p = CCodePrinter()
        return p.indent_code(code_lines)

    def update_boundary_conditions(self, code_dictionary):
        """ Generate OPSC code to affect a boundary condition update.

        :arg dict code_dictionary: The dictionary of OPSC code lines, with each key-value pair representing the code (value) for a particular stage (key) of the overall algorithm.
        :returns: The updated/modified code dictionary.
        :rtype: dict
        """
        from .bcs import ExchangeSelf
        from .kernel import Kernel

        bc_call = [[] for block in range(self.nblocks)]
        bc_exchange_code = [[] for block in range(self.nblocks)]

        for block in range(self.nblocks):
            for computation in self.boundary_condition[block].computations:
                if isinstance(computation, Kernel):
                    bc_call[block] += self.kernel_call(computation)
                elif isinstance(computation, ExchangeSelf):
                    call, code = self.bc_exchange_call_code(computation)
                    bc_call[block] += call
                    bc_exchange_code[block] += code
                else:
                    raise ValueError("Boundary condition of type %s cannot be classified" % (type(computation)))

        # Update the code dictionary
        code_dictionary['bc_exchange'] = '\n'.join(['\n'.join(bc_exchange_code[block]) for block in range(self.nblocks)])
        code_dictionary['bc_calls'] = '\n'.join(['\n'.join(bc_call[block]) for block in range(self.nblocks)])
        return code_dictionary

    def get_io(self, code_dictionary):
        """ As of now IO is performed only at the end of the simulation. No intermediate dumps are allowed.

        :arg dict code_dictionary: The dictionary of OPSC code lines, with each key-value pair representing the code (value) for a particular stage (key) of the overall algorithm.
        :returns: The updated/modified code dictionary.
        :rtype: dict
        """

        io_calls = [[] for block in range(self.nblocks)]
        io_time = [[] for block in range(self.nblocks)]
        for block in range(self.nblocks):
            # Process FileIO
            save_at = self.IO[block].save_after
            if len(save_at) == 1 and save_at[0] is True:
                name = self.simulation_parameters["name"] + '_' + str(int(self.simulation_parameters["niter"])) + '.h5'
                name = '\"' + name + '\"'
                io_calls[block] += self.hdf5_io(self.IO[block], name)
            else:
                name = self.simulation_parameters["name"] + '_' + str(int(self.simulation_parameters["niter"])) + '.h5'
                name = '\"' + name + '\"'
                io_calls[block] += self.hdf5_io(self.IO[block], name)
                # Time IO save at
                l = ccode(Mod('iteration+1', save_at[0]))
                calls = ['if(%s == 0)' % l] + [self.left_brace]
                # Character buffer array and name of the output
                calls += ['char buf[100];'] + ['sprintf(buf,\"%s_%%d.h5\",iteration);' % self.simulation_parameters["name"]]
                name = 'buf'
                calls += self.hdf5_io(self.IO[block], name) + [self.right_brace]
                io_time[block] += calls

        code_dictionary['io_calls'] = '\n'.join(['\n'.join(io_calls[block]) for block in range(self.nblocks)])
        code_dictionary['io_time'] = '\n'.join(['\n'.join(io_time[block]) for block in range(self.nblocks)])
        return code_dictionary

    def get_block_computation_kernels(self, instances):
        """ Get all computational kernel calls for each block.

        :returns: A list of lists, with each sublist containing all kernel calls for a particular block.
        :rtype: list of lists
        """
        # First process the inner time calls
        calls = [[] for block in range(self.nblocks)]
        for block in range(self.nblocks):
            for instance in instances[block]:
                if instance:
                    calls[block] += self.kernel_call(instance)
        return calls

    def kernel_call(self, computation):
        """ Generate an OPS kernel call via the ops_par_loop function.

        :arg computation: The computation to perform over the grid points.
        """
        iteration_range = self.iteration_range_name % self.iteration_range_index
        self.iteration_range_index += 1
        stencils = self.get_stencils(computation)
        kernel_calls = []

        # Iteration range (i.e. the range of the loops over the grid points).
        range_main = self.array('int', iteration_range, [r for ran in computation.ranges for r in ran])

        kernel_calls += ['ops_par_loop(%s, \"%s\", %s, %s, %s' % (computation.name, computation.computation_type, self.block_name, self.ndim, iteration_range)]

        # Do the grid-based inputs first.
        grid_based = [self.ops_argument_call(inp, stencils[inp], self.dtype, self.ops_access['inputs'])
                      for inp in computation.inputs.keys() if inp.is_grid] + \
                     [self.ops_argument_call(inp, stencils[inp], self.dtype, self.ops_access['outputs'])
                      for inp in computation.outputs.keys() if inp.is_grid] + \
                     [self.ops_argument_call(inp, stencils[inp], self.dtype, self.ops_access['inputoutput'])
                      for inp in computation.inputoutput.keys() if inp.is_grid]
        # Globals
        nongrid = [self.ops_global_call(inp, value, self.dtype, self.ops_access['inputs'])
                   for inp, value in computation.inputs.iteritems() if not inp.is_grid] + \
                  [self.ops_global_call(inp, value, self.dtype, self.ops_access['outputs'])
                   for inp, value in computation.outputs.iteritems() if not inp.is_grid] + \
                  [self.ops_global_call(inp, value, self.dtype, self.ops_access['inputoutput'])
                   for inp, value in computation.inputoutput.iteritems() if not inp.is_grid]
        # Reductions
        if computation.reductions:
            nongrid += [self.ops_argument_reduction(inp, self.ops_access['reduction'])
                        for inp in computation.reductions if isinstance(inp, ReductionVariable)]

        if computation.has_Idx:
            nongrid.append(self.grid_index_call())

        kernel_calls = kernel_calls + grid_based + nongrid
        call = [k + ',' for k in kernel_calls[:-1]]
        call = [range_main] + call + [kernel_calls[-1] + self.right_parenthesis + self.end_of_statement] + ['\n']

        return call

    def get_stencils(self, computation):
        stencils = {}
        dicts = [computation.inputs, computation.outputs, computation.inputoutput]
        for d in dicts:
            for key, value in d.iteritems():
                if key.is_grid:
                    relative_stencil = self.relative_stencil(value)
                    stencil = self.list_to_string(relative_stencil)
                    if stencil not in self.stencil_dictionary.keys():
                        self.stencil_dictionary[stencil] = self.stencil_name % self.stencil_number
                        self.stencil_number += 1
                    # Update the stencils to be returned.
                    stencils[key] = self.stencil_dictionary[stencil]
        return stencils

    def list_to_string(self, l):
        """ Convert a list to a string, where each list element is separated by a comma.

        :arg list l: The list to convert into a string.
        :returns: The list elements converted to a string, where each list element is separated by a comma.
        :rtype: str
        """

        s = ','.join([str(element) for element in l])
        return s

    def relative_stencil(self, value):
        """ Returns the relative stencil wrt the grid location. i.e. grid indices eg(i0,i1,i2) are replaced with (0,0,0). """

        if isinstance(value, list):
            pass
        else:
            value = [value]

        indices = []
        for val in value:
            indices.append(tuple([ind for ind in val if ind != EinsteinTerm('t')]))

        return_val = []
        for va in indices:
            out = []
            for number, v in enumerate(va):
                outv = v
                for a in v.atoms(Symbol):
                    outv = outv.subs(a, 0)
                out.append(outv)
            return_val.append(out)

        return_val = self.sort_stencil_indices(return_val)

        return return_val

    def sort_stencil_indices(self, indexes):
        """ Helper function for relative_stencil. Sorts the relative stencil. """
        if len(indexes[0]) > 1:
            for dim in range(len(indexes[0])):
                indexes = sorted(indexes, key=lambda indexes: indexes[dim])
            temp = flatten(list(list(t) for t in indexes))
        else:
            # print(indexes)
            indexes = [sorted(indexes)]
            temp = flatten(list(t) for t in indexes)
        return temp

    def grid_index_call(self):
        """ The call to the OPS helper function to get the grid point index.

        :returns: The call to ops_arg_idx().
        :rtype: str
        """
        return 'ops_arg_idx()'

    def ops_global_call(self, array, indices, precision, access_type):
        arr = array[tuple(indices[0])]
        template = 'ops_arg_gbl(&%s, %d, \"%s\", %s)'
        return template % (arr, 1, self.dtype, access_type)

    def ops_argument_call(self, array, stencil, precision, access_type):
        template = 'ops_arg_dat(%s, %d, %s, \"%s\", %s)'
        return template % (array, 1, stencil, self.dtype, access_type)

    def ops_argument_reduction(self, name, access_type):
        template = 'ops_arg_reduce(%s, %d, \"%s\", %s)'
        return template % (name, 1, self.dtype, access_type)

    def bc_exchange_call_code(self, instance):
        off = 0
        halo = 'halo'
        # Name of the halo exchange
        name = self.halo_exchange_name % (self.halo_exchange_number)
        self.halo_exchange_number = self.halo_exchange_number + 1
        code = ['%s Boundary condition exchange code' % self.line_comment]
        code += ['ops_halo_group %s %s' % (name, self.end_of_statement)]
        code += [self.left_brace]
        code += ['int halo_iter[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_size]), self.end_of_statement)]
        code += ['int from_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_from]), self.end_of_statement)]
        code += ['int to_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_to]), self.end_of_statement)]
        # dir in OPSC. FIXME: Not sure what it is, but 1 to ndim works.
        code += ['int dir[] = {%s}%s' % (', '.join([str(ind+1) for ind in range(len(instance.transfer_to))]), self.end_of_statement)]
        # Process the arrays
        for arr in instance.transfer_arrays:
            code += ['ops_halo %s%d = ops_decl_halo(%s, %s, halo_iter, from_base, to_base, dir, dir)%s'
                     % (halo, off, arr.base, arr.base, self.end_of_statement)]
            off = off+1
        code += ['ops_halo grp[] = {%s}%s' % (','.join([str('%s%s' % (halo, of)) for of in range(off)]), self.end_of_statement)]
        code += ['%s = ops_decl_halo_group(%d,grp)%s' % (name, off, self.end_of_statement)]
        code += [self.right_brace]
        # Finished OPS halo exchange, now get the call
        call = ['%s Boundary condition exchange calls' % self.line_comment, 'ops_halo_transfer(%s)%s' % (name, self.end_of_statement)]
        return call, code

    def initialise_dat(self):
        """ Initialise OPS dats (i.e. datasets).

        :returns: The declaration code in OPSC format.
        :rtype: str
        """

        code = ['%s Initialise/allocate OPS dataset.' % (self.line_comment)]
        dtype_int = 'int'
        if not self.multiblock:
            grid = self.grid[0]
            code += [self.array(dtype_int, 'halo_p', [halo[1] for halo in grid.halos])]
            code += [self.array(dtype_int, 'halo_m', [halo[0] for halo in grid.halos])]
            code += [self.array(dtype_int, 'size', grid.shape)]
            code += [self.array(dtype_int, 'base', [0 for g in grid.shape])]
            code += ['%s* val = NULL;' % (self.dtype)]
            init_format = '%%s = ops_decl_dat(%s, 1, size, base, halo_m, halo_p, val, \"%%s\", \"%%s\")%s' % (self.block_name, self.end_of_statement)
            inits = [init_format % (arr, self.dtype, arr) for arr in self.grid_based_arrays]
            code = code + inits
        else:
            raise NotImplementedError("Multi-block is not implemented")
        return code

    def declare_stencils(self):
        """ Declare all the stencils used in the code. We do not differentiate between the stencils for each block.

        :returns: The OPSC code declaring the stencil.
        :rtype: str
        """

        code = ['%s Declare all the stencils used ' % (self.line_comment)]
        dtype_int = 'int'
        sten_format = 'ops_stencil %%s = ops_decl_stencil(%%d,%%d,%%s,\"%%s\")%s' % (self.end_of_statement)
        for key, value in self.stencil_dictionary.iteritems():
            count = len(key.split(',')) / self.ndim
            # 'value' is the name in the stencil's format
            code += [self.array(dtype_int, value + "_temp", [key])]
            code += [sten_format % (value, self.ndim, count, value + "_temp", key)]
        return code

    def hdf5_io(self, instance, name):
        code = []
        block_to_hdf5 = ["ops_fetch_block_hdf5_file(%s, %s)%s" % (self.block_name, name, self.end_of_statement)]
        code += block_to_hdf5
        # Then write out each field.
        for c in instance.save_arrays:
            variables_to_hdf5 = ["ops_fetch_dat_hdf5_file(%s, %s)%s" % (c, name, self.end_of_statement)]
            code += variables_to_hdf5
        return code

    def get_block_computations(self):
        """ Get all the block computations to be performed.
        Extra stuff like diagnostic computations or boundary condition computations should be added here.

        :returns: A list of computations written in OPSC format.
        :rtype: list
        """
        from .kernel import Kernel

        kernels = [[] for block in range(self.nblocks)]
        for block in range(self.nblocks):
            # Get all the computations to be performed. Add computations as needed.
            block_computations = []
            if self.spatial_discretisation[block].computations:
                block_computations += self.spatial_discretisation[block].computations
            if self.temporal_discretisation[block].computations:
                block_computations += self.temporal_discretisation[block].computations
            if self.temporal_discretisation[block].start_computations:
                block_computations += self.temporal_discretisation[block].start_computations
            if self.temporal_discretisation[block].end_computations:
                block_computations += self.temporal_discretisation[block].end_computations
            if self.initial_conditions[block].computations:
                block_computations += self.initial_conditions[block].computations
            if self.diagnostics:
                for inst in self.diagnostics[block]:
                    block_computations += inst.computations
            if self.boundary_condition[block].computations:
                block_computations += [t for t in self.boundary_condition[block].computations if isinstance(t, Kernel)]

            for computation in block_computations:
                kernels[block] += self.kernel_computation(computation, block)
        return kernels

    def kernel_computation(self, computation, block_number):
        """ Generate the kernel for the computation. This acts as a helper function for the block computations.

        :arg computation: The computation to perform over the grid points in a particular block.
        :arg block_number: The number of the block over which to perform the computation.
        :returns: A list of OPSC code lines performing the computation.
        :rtype: list
        """
        from .grid import GridVariable

        header = []

        if computation.name is None:
            computation.name = self.computational_kernel_names[block_number] % self.kernel_name_number[block_number]

        # Indexed objects based on the grid that are inputs/outputs or inouts. This is used to write the pointers to the kernel.
        # Grid-based objects
        grid_based = ([self.ops_header['inputs'] % (self.dtype, inp) for inp in computation.inputs.keys() if inp.is_grid] +
                      [self.ops_header['outputs'] % (self.dtype, inp) for inp in computation.outputs.keys() if inp.is_grid] +
                      [self.ops_header['inputoutput'] % (self.dtype, inp) for inp in computation.inputoutput.keys() if inp.is_grid])

        # Non grid-based objects
        nongrid = ([self.ops_header['inputs'] % (self.dtype, inp) for inp in computation.inputs.keys() if not inp.is_grid] +
                   [self.ops_header['outputs'] % (self.dtype, inp) for inp in computation.outputs.keys() if not inp.is_grid] +
                   [self.ops_header['inputoutput'] % (self.dtype, inp) for inp in computation.inputoutput.keys() if not inp.is_grid])

        header += grid_based + nongrid

        if computation.reductions:
            header += [self.ops_header['reduction'] % (self.dtype, inp) for inp in computation.reductions]
        if computation.has_Idx:
            header += [self.ops_header['Idx'] % ('idx')]

        header = ['void ' + computation.name + self.left_parenthesis + ' , '.join(header) + self.right_parenthesis]
        header += [self.left_brace]
        code = header
        ops_accs = self.get_OPS_ACC_number(computation)

        for equation in computation.equations:
            code_kernel, self.rational_constants = ccode(equation, ops_accs, self.rational_constants)
            if isinstance(equation.lhs, GridVariable):

                code += [self.dtype + ' ' + code_kernel + self.end_of_statement]
            else:
                code += [code_kernel + self.end_of_statement]

        code += [self.right_brace] + ['\n']

        self.update_definitions(computation)

        # Update the kernel name index
        self.kernel_name_number[block_number] += 1

        return code

    def write_computational_routines(self, kernels):
        """ Write the computational routines to files. """
        for block in range(self.nblocks):
            # Write the code for each block.
            code_lines = ["#ifndef block_%d_KERNEL_H" % block + '\n' + "#define block_%d_KERNEL_H" % block + '\n']
            code_lines += kernels[block]
            code_lines += ["#endif"]
            kernel_file = open(self.CODE_DIR + '/' + self.computational_routines_filename[block], 'w')
            kernel_file.write('\n'.join(code_lines))
            kernel_file.close()
        return

    def loop_open(self, var, range_of_loop):
        """ The head of a for-loop.

        :arg str var: The name of the iteration variable.
        :arg range_of_loop: A list or tuple containing the start and end points of the loop over 'var'.
        :returns: The head of the loop in OPSC format.
        :rtype: str
        """
        return 'for (int %s=%d; %s<%d; %s++)%s' % (var, range_of_loop[0], var, range_of_loop[1], var, self.left_brace)

    def loop_close(self):
        """ Close a for loop.

        :returns: A closing brace to close a for loop.
        :rtype: str
        """
        return self.right_brace

    def header(self):
        """ Header code.

        :returns: A list of header lines in OPSC format.
        :rtype: list
        """

        code = []
        code += ['#include <stdlib.h>']
        code += ['#include <string.h>']
        code += ['#include <math.h>']
        code += ['%s Global constants in the equations are' % self.line_comment]
        for constant in self.constants:
            if isinstance(constant, IndexedBase):
                code += ['%s %s[%d]%s' % (self.dtype, constant, constant.ranges, self.end_of_statement)]
            elif isinstance(constant, str):
                code += ['%s %s%s' % (self.dtype, constant, self.end_of_statement)]
            elif constant.is_integer:
                code += ['int %s%s' % (constant, self.end_of_statement)]
            else:
                code += ['%s %s%s' % (self.dtype, constant, self.end_of_statement)]

        # Include constant declaration
        code += ['// OPS header file']
        code += ['#define OPS_%sD' % self.ndim]
        code += ['#include "ops_seq.h"']
        # Include the kernel file names
        code += ['#include "%s"' % name for name in self.computational_routines_filename]
        return code

    def main_start(self):
        """ The famous 'int main' statement. """

        return ['%s main program start' % self.line_comment, 'int main (int argc, char **argv) ', self.left_brace]

    def ops_init(self, diagnostics_level=None):
        """ The default diagnostics level is 1, which offers no diagnostic information and should be used for production runs.
        Refer to OPS user manual for more information.

        :arg int diagnostics_level: The diagnostics level. If None, the diagnostic level defaults to 1.
        :returns: The call to ops_init.
        :rtype: list
        """
        out = ['%s Initializing OPS ' % self.line_comment]
        if diagnostics_level:
            self.ops_diagnostics = True
            return out + ['ops_init(argc,argv,%d)%s' % (diagnostics_level, self.end_of_statement)]
        else:
            self.ops_diagnostics = False
            return out + ['ops_init(argc,argv,%d)%s' % (1, self.end_of_statement)]

    def ops_diagnostics(self):
        # WARNING: Untested OPS diagnostics output. Need to check if this gives the correct result or not.
        if self.ops_diagnostics:
            return ['ops diagnostic output()']
        else:
            return []
        return

    def ops_partition(self):
        """ Initialise an OPS partition for the purpose of multi-block and/or MPI partitioning.

        :returns: The partitioning code in OPSC format. Each line is a separate list element.
        :rtype: list
        """

        return ['%s Init OPS partition' % self.line_comment, 'ops_partition(\"\")%s' % self.end_of_statement]

    def ops_timers(self):
        """ OPS timer variable declaration. """
        start = ["cpu_start", "elapsed_start"]
        end = ["cpu_end", "elapsed_end"]
        timer_start = ["double %s, %s%s" % (start[0], start[1], self.end_of_statement)] + ["ops_timers(&%s, &%s)%s" % (start[0], start[1], self.end_of_statement)]
        timer_end = ["double %s, %s%s" % (end[0], end[1], self.end_of_statement)] + ["ops_timers(&%s, &%s)%s" % (end[0], end[1], self.end_of_statement)]
        timing_eval = self.ops_print_timings(start, end)
        return timer_start, timer_end, timing_eval

    def ops_print_timings(self, start, end):
        """ Generate OPSC code to print out run-time information.

        :returns: A list of code lines for printing out run-time information generated by OPS timers.
        :rtype: list
        """

        code = []
        code += ["ops_printf(\"\\nTimings are:\\n\")%s" % self.end_of_statement]
        code += ["ops_printf(\"-----------------------------------------\\n\")%s" % self.end_of_statement]
        code += ["ops_printf(\"Total Wall time %%lf\\n\",%s-%s)%s" % (end[1], start[1], self.end_of_statement)]
        return code

    def define_block(self):
        """ Define an OPS block. """
        code = ['%s Defining block in OPS Format' % (self.line_comment)]
        if not self.multiblock:
            # No dynamic memory allocation required
            code += ['ops_block %s%s' % (self.block_name, self.end_of_statement)]
        else:
            code += ['ops_block *%s = (ops_block *)malloc(%s*sizeof(ops_block*))%s' % (self.block_name, self.nblocks, self.end_of_statement)]
        # print('\n'.join(code))
        return code

    def initialise_block(self):
        """ Initialise an OPS block. """
        code = ['%s Initialising block in OPS Format' % (self.line_comment)]
        if not self.multiblock:
            code += ['%s = ops_decl_block(%d, \"%s\")%s' % (self.block_name, self.ndim, self.block_name, self.end_of_statement)]
        else:
            raise NotImplementedError("Multi-block is not implemented")
        # print('\n'.join(code))
        return code

    def define_dat(self):
        """ Define OPS dats (i.e. datasets) for all grid-based arrays. """

        code = ['%s Define dataset' % (self.line_comment)]
        if not self.multiblock:
            def_format = 'ops_dat %%s%s' % self.end_of_statement
            code += [def_format % arr for arr in self.grid_based_arrays]
        else:
            raise NotImplementedError("Multi-block is not implemented")
        return code

    def ops_exit(self):
        """ Calls ops_exit. """
        return ['%s Exit OPS ' % self.line_comment, 'ops_exit()%s' % self.end_of_statement]

    def footer(self):
        """ Footer code in OPSC format. """
        code = self.ops_exit()
        return code

    def check_consistency(self, grid, spatial_discretisation, temporal_discretisation, boundary_condition, initial_conditions, IO):
        """ Check the consistency of the inputs. """

        self.grid = self.to_list(grid)
        length = len(self.grid)

        self.spatial_discretisation = self.to_list(spatial_discretisation)
        if len(self.spatial_discretisation) != length:
            raise AlgorithmError("The length of spatial solution does not match the grid.")

        self.temporal_discretisation = self.to_list(temporal_discretisation)
        if len(self.temporal_discretisation) != length:
            raise AlgorithmError("The length of temporal solution does not match the grid.")

        self.boundary_condition = self.to_list(boundary_condition)
        if len(self.boundary_condition) != length:
            raise AlgorithmError("The length of boundary does not match the grid.")

        self.initial_conditions = self.to_list(initial_conditions)
        if len(self.initial_conditions) != length:
            raise AlgorithmError("The length of initial_conditions does not match the grid.")

        self.IO = self.to_list(IO)
        if len(self.IO) != length:
            raise AlgorithmError("The length of IO does not match the grid.")

        return

    def to_list(self, var):
        """ Convert a non list object into a list.

        :arg var: The non list variable. Note: if this is a list, then it is simply returned immediately without modification.
        :returns: The non-list variable converted into a list.
        :rtype: list
        """
        if isinstance(var, list):
            return var
        else:
            return [var]

    def update_definitions(self, computation):
        """ Update the grid based arrays and constants to be declared. """

        arrays = set([inp for inp in computation.inputs.keys() if inp.is_grid] +
                     [inp for inp in computation.outputs.keys() if inp.is_grid] +
                     [inp for inp in computation.inputoutput.keys() if inp.is_grid])
        constant_arrays = set([inp for inp in computation.inputs.keys() if not inp.is_grid] +
                              [inp for inp in computation.outputs.keys() if not inp.is_grid] +
                              [inp for inp in computation.inputoutput.keys() if not inp.is_grid])
        constants = set(computation.constants)

        # add the symbols in the range of evaluation to constants
        const = set(flatten([list(r.atoms(Symbol)) for ro in computation.ranges for r in ro
                             if not isinstance(r, int)]))

        self.grid_based_arrays = self.grid_based_arrays.union(arrays)

        self.constants = self.constants.union(constant_arrays).union(constants).union(self.rational_constants.values())\
            .union(const)
        # update the simulation parameters as these are used for writing the constants
        self.simulation_parameters.update(dict(zip([str(v) for v in self.rational_constants.values()],
                                                   list(self.rational_constants.keys()))))

        # Update reduction variables
        self.reduction_variables = self.reduction_variables.union(set(computation.reductions))
        return

    def get_OPS_ACC_number(self, computation):
        """ Helper function for writing OPS kernels, which obtains all the of the OPS_ACCs.

        :arg computation: The computational kernel to write.
        :returns: A dictionary of OPS_ACC's. """

        ops_accs = {}
        allidbs = list(computation.inputs.keys()) + list(computation.outputs.keys()) + list(computation.inputoutput.keys())
        grid_based = [al for al in allidbs if al.is_grid]
        # All grid-based OPS_ACCs
        for no, inp in enumerate(grid_based):
            ops_accs[inp] = 'OPS_ACC%d' % no
        # Non grid-based stuff
        nongrid = set(allidbs).difference(set(grid_based))
        for no, inp in enumerate(nongrid):
            ops_accs[inp] = None
        return ops_accs

    def array(self, dtype, name, values):
        """ Declare inline arrays in OPSC/C

        :arg dtype: The data type of the array.
        :arg name: The name of the array.
        :arg size: The size of the array.
        :arg vals: The list of values.
        :returns: The
        :rtype: str
        """
        return '%s %s[] = {%s}%s' % (dtype, name, ', '.join([str(s) for s in values]), self.end_of_statement)

    def translate(self):
        # Translate the generated code using the OPSC translator.
        LOG.debug("Translating OPSC code...")
        exit_code = subprocess.call("cd %s; python -c 'from ops_translator.c import ops as translator; translator.main([\"%s.cpp\"])'" % (self.CODE_DIR, self.simulation_parameters["name"]), shell=True)
        if(exit_code != 0):
            # Something went wrong
            LOG.error("Unable to translate OPSC code. Check that OPS is installed.")
        return


class AlgorithmError(Exception):

    """ An Exception that occurs  """

    pass

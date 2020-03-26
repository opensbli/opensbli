"""@brief Simulation monitors to report data values during the simulation.
   @author David J. Lusher
   @details Can be used to generate a large number of flow samples to construct a time signal.
"""
from sympy import flatten

class Monitor(object):
    def __init__(self, flow_var, probe_loc, numbering):
        self.flow_var = flow_var
        self.probe_loc = probe_loc
        self.probe_no = numbering
        return

class SimulationMonitor(object):
    def __init__(self, arrays, probe_locations, block, print_frequency = 250, fp_precision=10, NaNcheck=True, output_file=None):
        """ Class to enable access of dataset values during the simulation.
        :arg list arrays: A list of DataSets to monitor during the simulation.
        :arg list probe_locations: A list of tuples giving the (i,j,k) grid index location of the probe.
        :arg object block: An OpenSBLI simulation block.
        :arg int print_frequency: The iteration frequency at which to print the output.
        :arg int fp_precision: The number of decimal places to format the output.
        :arg bool NaNcheck: Adds an optional call to ops_NaNcheck().
        :arg str output_file: The option to write the output directly to a log file, otherwise defaults to stdout."""

        # Check number of probes equals the number of input arrays
        if len(arrays) != len(probe_locations):
            raise ValueError("The number of arrays must equal the number of probe locations.")
        self.monitors = [Monitor(var, loc, index) for index, (var,loc) in enumerate(zip(arrays, probe_locations))]
        self.components = []
        self.frequency = print_frequency
        self.fp_precision = fp_precision
        self.output_file = output_file
        self.block = block
        self.ndim = block.ndim
        self.filename = 'reductions.h'
        self.NaNcheck = NaNcheck
        return

    def add_components(self, components):
        """ Adds the given components to the simulation monitoring
        :param components: the components to be added to the timers, this can be a list or an individual component
        :return: None. """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]
        return

    def write_latex(self, latex):
        latex.write_string("Simulation monitoring start\\\\\n")
        for c in self.components:
            c.write_latex(latex)
        latex.write_string("Simulation monitoring end\\\\\n")
        return

    @property
    def opsc_code(self):
        code = self.opsc_start
        code += self.opsc_middle
        code += self.opsc_end
        return code

    def initial_declarations(self, M):
        """ Initial declarations for the OPS data access."""
        name, number = str(M.flow_var), M.probe_no
        declarations = ["// Monitoring of %s" % name]
        declarations += ["ops_reduction reduce_%d_%s = ops_decl_reduction_handle(sizeof(double), \"double\", \"reduction_%d_%s\");" % (number, name, number, name)]
        declarations += ["double %s_%d_output = 0.0;" % (name, number)]
        return declarations

    def reduction_range(self, M, idx):
        """ Sets up the grid locations and iteration range."""
        name, number = str(M.flow_var), M.probe_no
        if self.ndim == 1:
            i = M.probe_loc
            return ["int i%d = %s;" % (idx, str(i))] + ["int monitor_range_%d_%s[] = {i%d, i%d+1};" % (number, name, idx, idx)]
        elif self.ndim == 2:
            i, j = M.probe_loc
            return ["int i%d = %s, j%d = %s;" % (idx, str(i), idx, str(j))] + ["int monitor_range_%d_%s[] = {i%d, i%d+1, j%d, j%d+1};" % (number, name, idx, idx, idx, idx)]
        elif self.ndim == 3:
            i, j, k = M.probe_loc
            return ["int i%d = %s, j%d = %s, k%d = %s;" % (idx, str(i), idx, str(j), idx, str(k))] + ["int monitor_range_%d_%s[] = {i%d, i%d+1, j%d, j%d+1, k%d, k%d+1};" % (number, name, idx, idx, idx, idx, idx, idx)]

    def par_loop_declaration(self, M):
        """ Defines the parallel loop templates for the reductions."""
        name, number = str(M.flow_var), M.probe_no
        output_code = ["ops_par_loop(monitor_%d_%s, \"Reduction %s_%d\", %s, %d, monitor_range_%d_%s," % (number, name, name, number, self.block.blockname, self.ndim, number, name)]
        output_code += ["ops_arg_dat(%s, 1, zero_stencil, \"double\", OPS_READ)," % name]
        output_code += ["ops_arg_reduce(reduce_%d_%s, 1, \"double\", OPS_INC));" % (number, name)]
        return output_code

    def generate_result(self, M):
        """ Obtains the result of the OPS reduction."""
        name, number = str(M.flow_var), M.probe_no
        return ["ops_reduction_result(reduce_%d_%s, &%s_%d_output);\n" % (number, name, name, number)]

    def generate_kernel_code(self, M):
        """ Generates the kernel definitions to be added to the reductions header file."""
        name, number = str(M.flow_var), M.probe_no
        code = ["void monitor_%d_%s(const double *%s, double *reduce_%d_%s){\n" % (number, name, name, number, name)]
        indices = ','.join(['0' for _ in range(self.ndim)])
        code += ["*reduce_%d_%s = %s[OPS_ACC0(%s)];\n}" % (number, name, name, indices)] + ["\n\n"]
        return code

    @property
    def add_NaN_check(self):
        return ['ops_NaNcheck(%s);' % str(self.monitors[0].flow_var)]

    @property
    def write_reductions_file(self):
        """ Creates a new header file to define the reduction kernels."""
        code = []
        f = open(self.filename, 'w')
        f.write("#ifndef REDUCTIONS_H\n")
        f.write("#define REDUCTIONS_H\n")
        for M in self.monitors:
            f.write(''.join(self.generate_kernel_code(M)))
        f.write('#endif\n')
        f.close()
        return

    @property
    def format_output(self):
        """ Controls the printing format for the output."""
        placeholders = ', '.join(["%d"] + ["%%.%df" % self.fp_precision for _ in range(len(self.monitors)+1)])
        iterations = ['iter+1', '(iter+1)*dt']
        variables = ["%s_%d_output" % (str(M.flow_var), M.probe_no) for M in self.monitors]
        # Normalise mean quantities by the number of iterations
        for i, var in enumerate(variables):
            if 'mean' in var:
                variables[i] = '%s/%s' % (var, '(iter+1)')
        variables = ', '.join(iterations + variables)
        # Check if the output should be written directly to a log file
        if self.output_file:
            output_print = ["ops_fprintf(f, \"%s\\n\", %s);" % (placeholders, variables)]
        else:
            output_print = ["ops_printf(\"%s\\n\", %s);" % (placeholders, variables)]
        return ["// Write the output values"] + output_print

    @property
    def generate_reduction_loops(self):
        """ Creates a block of code for each flow variable being monitored."""
        output_code = []
        for index, M in enumerate(self.monitors):
            output_code += self.initial_declarations(M)
            output_code += self.reduction_range(M, index)
            output_code += self.par_loop_declaration(M)
            output_code += self.generate_result(M)
        return output_code

    @property
    def initial_print(self):
        """ Prints the headers at the top of the output once at the start."""
        headers = ['Iteration', 'Time'] + ['%s%s' % (str(M.flow_var), str(M.probe_loc)) for M in self.monitors]
        headers = ', '.join(headers)
        if self.output_file:
            output_code = ["if (iter == 0){\nops_fprintf(f, \"%s\\n\");}" % (headers)]
        else:
            output_code = ["if (iter == 0){\nops_printf(\"%s\\n\");}" % (headers)]
        return output_code

    @property
    def declare_stencils(self):
        """ Zero stencil declaration as no relative indexing is being used at present."""
        if self.ndim == 1:
            return ["int zero_stencil_indices[] = {0};"] + ["ops_stencil zero_stencil = ops_decl_stencil(1, 1, zero_stencil_indices, \"zero_stencil\");\n"]
        elif self.ndim == 2:
            return ["int zero_stencil_indices[] = {0,0};"] + ["ops_stencil zero_stencil = ops_decl_stencil(2, 1, zero_stencil_indices, \"zero_stencil\");\n"]
        elif self.ndim == 3:
            return ["int zero_stencil_indices[] = {0,0,0};"] + ["ops_stencil zero_stencil = ops_decl_stencil(3, 1, zero_stencil_indices, \"zero_stencil\");\n"]
    
    @property
    def opsc_start(self):
        starting_code = ["// Data access for simulation monitoring"]
        if self.NaNcheck:
            starting_code += self.add_NaN_check
        starting_code += self.initial_print + self.declare_stencils 
        return starting_code

    @property
    def opsc_middle(self):
        middle_code = []
        middle_code += self.generate_reduction_loops
        return middle_code
    
    @property
    def opsc_end(self):
        end_code = self.format_output
        self.write_reductions_file
        return end_code

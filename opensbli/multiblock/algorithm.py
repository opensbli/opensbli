from opensbli.multiblock.blockcollection import MultiBlock as MB

from sympy import flatten, pprint, Idx, Equality
from opensbli.code_generation.latex import LatexWriter
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.equation_types.opensbliequations import SimulationEquations, NonSimulationEquations, ConstituentRelations
from opensbli.equation_types.metric import MetricsEquation
from opensbli.core.opensbliobjects import Constant, DataSetBase, ConstantObject
from opensbli.core.datatypes import Int
from opensbli.code_generation.algorithm.common import BeforeSimulationStarts, AfterSimulationEnds, InTheSimulation
import copy


class Loop(object):
    """ Base object representing loops in an algorithm
    """
    pass


class MainPrg(Loop):
    """ Main program of the algorithm. The attribute components are used to loop over
    """

    def __init__(self):
        self.components = []
        return

    def __str__(self):
        return "%s" % (self.__class__.__name__)

    def add_components(self, components):
        """ Adds the given components to the main program

        :param components: the components to be added to the main program, this can be a list or an individual component
        :return: None
        """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]
        return

    def write_latex(self, latex):
        """ Writes the LaTeX of the main program by looping over the components

        :param latex: the opened file pointer to which the generated LaTeX code should be written.
        :returns: None

        """
        latex.write_string("Starting of the main program\\\\\n")
        for c in self.components:
            c.write_latex(latex)
        latex.write_string("End of Main")
        return

    @property
    def opsc_code(self):
        """ Writes the OPS C version of the code by looping over the components
        """
        code = []
        code += [self.opsc_start]
        for c in flatten(self.components):
            code += c.opsc_code
        code += [self.opsc_end]
        return code

    @property
    def opsc_start(self):
        """ Starting the program in OPS C
        """
        return "int main(int argc, char **argv) \n{"

    @property
    def opsc_end(self):
        """ Ending the program in OPS C
        """
        return "//Main program end \n}"


class Condition(object):
    """ Used for setting the conditions in the program, example
    """

    def __init__(self, condition):
        self.condition = condition
        self.components = []
        return

    def add_components(self, components):
        """ Adds the given components to the Condition

        :param components: the components to be added to the Condition, this can be a list or an individual component
        :return: None
        """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]

    def write_latex(self, latex):
        """ Writes the LaTeX of the Condition by looping over the components

        :param latex: the opened file pointer to which the generated LaTeX code should be written.
        :returns: None

        """
        latex.write_string("Condition %s" % (self.condition))
        for c in self.components:
            c.write_latex(latex)
        latex.write_string("Ending condition loop %s")
        return

    @property
    def opsc_code(self):
        """ Writes the OPS C version of the code by looping over the components
        """
        code = self.opsc_condition_start
        for c in self.components:
            code += c.opsc_code
        code += self.opsc_condition_end
        return code

    @property
    def opsc_condition_start(self):
        """ The starting loop for a if condition in OPS C
        """
        from opensbli.code_generation.opsc import OPSCCodePrinter
        code = OPSCCodePrinter().doprint(self.condition)
        code = 'if (' + code + '){'
        return [code]

    @property
    def opsc_condition_end(self):
        """ The loop end for an if condition in OPS C
        """
        return ['}']


class DoLoop(Loop):
    def __init__(self, iterator):
        self.loop = iterator
        self.components = []
        return

    def add_components(self, components):
        """ Adds the given components to the Do Loop

        :param components: the components to be added to the DoLoop, this can be a list or an individual component
        :return: None
        """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]
        return

    def write_latex(self, latex):
        """ Writes the LaTeX of the components of the Do loop

        :param latex: the opened file pointer to which the generated LaTeX code should be written.
        :returns: None

        """
        latex.write_string("Do loop %s,%s, %s\\\\\n" % (self.loop, self.loop.lower, self.loop.upper))
        for c in self.components:
            c.write_latex(latex)
        latex.write_string("Ending Do loop %s\\\\\n" % (self.loop))
        return

    @property
    def opsc_code(self):
        """ Writes the OPS  version of the Do loop with the components
        """
        code = []
        code += [self.opsc_start]
        for c in self.components:
            code += c.opsc_code
        code += [self.opsc_end]
        return code

    @property
    def opsc_start(self):
        """ Do loop in OPSC is a for loop, and the starting of the for loop is written by this function
        """
        return "for(int %s=%s; %s<=%s; %s++)\n{" % (self.loop, str(self.loop.lower), self.loop, str(self.loop.upper), self.loop)

    @property
    def opsc_end(self):
        """ Ending of a for loop in OPSC
        """
        return "}"


class DefDecs(object):
    """ Definitions and declarations in a program. This write latex and OPS C code printing functions are
    not used currently but these will be added in the future releases
    """

    def __init__(self):
        self.components = []

    def add_components(self, components):
        """ Adds the given components to the definitions and declarations

        :param components: the components to be added to the DefDecs, this can be a list or an individual component
        :return: None
        """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]

    def write_latex(self, latex):
        """ Writes the LaTeX of the definitions are declarations

        :param latex: the opened file pointer to which the generated LaTeX code should be written.
        :returns: None

        .. note:: This is not used currently

        """
        for c in self.components:
            if isinstance(c, Constant):
                if c.is_input:
                    l = "Is an input"
                else:
                    l = "Is not an input"
                latex.write_string("DefDec constant %s, %s\\\\\n" % (str(c), l))
            elif isinstance(c, DataSetBase):
                latex.write_string("DefDec dataset %s\\\\\n" % (str(c)))
        return

    @property
    def opsc_code(self):
        """ Writes the Definitions and declarations in OPSC,

        .. note:: This is not used currently
        """
        code = []
        for c in self.components:
            if isinstance(c, Constant):
                code += ["DefDec constant %s" % (str(c))]
            elif isinstance(c, DataSetBase):
                code += ["DefDec dataset %s" % (str(c))]
        return code


class Timers(object):
    def __init__(self, number):
        self.components = []
        self.number = number
        return

    @property
    def _start_variables(self):
        return ["cpu_start%d" % (self.number), "elapsed_start%d" % (self.number)]

    @property
    def _end_variables(self):
        return ["cpu_end%d" % (self.number), "elapsed_end%d" % (self.number)]

    def add_components(self, components):
        """ Adds the given components to the timers

        :param components: the components to be added to the timers, this can be a list or an individual component
        :return: None
        """
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]
        return

    def write_latex(self, latex):
        """ Writes the LaTeX of the times

        :param latex: the opened file pointer to which the generated LaTeX code should be written.
        :returns: None

        """
        latex.write_string("Initialising timers\\\\\n")
        for c in self.components:
            c.write_latex(latex)
        latex.write_string("Ending timers\\\\\n")
        return

    @property
    def opsc_code(self):
        """ OPSC code for the timers, this combines the codes from different components

        .. note:: The data type for timer variables `double`
        """
        code = []
        code += self.opsc_start_timer
        for c in self.components:
            code += c.opsc_code
        code += self.opsc_end_timer
        return code

    @property
    def opsc_start_timer(self):
        """ OPSC code for the starting the timers, this is called from Timers.opsc_code()
        """
        start = self._start_variables
        timer_start = ["double %s, %s;" % (start[0], start[1])] + ["ops_timers(&%s, &%s);" % (start[0], start[1])]
        return timer_start

    @property
    def opsc_end_timer(self):
        """ OPSC code for the ending the timers, this is called from Timers.opsc_code()
        """
        end = self._end_variables
        code = []
        code += ["double %s, %s;" % (end[0], end[1])] + ["ops_timers(&%s, &%s);" % (end[0], end[1])]
        code += self.timing_result_opsc
        return code

    @property
    def timing_result_opsc(self):
        """ Generates the code for printing the timing results in OPSC
        """
        code = []
        code += ["ops_printf(\"\\nTimings are:\\n\");"]
        code += ["ops_printf(\"-----------------------------------------\\n\");"]
        code += ["ops_printf(\"Total Wall time %%lf\\n\",%s-%s);" % (self._end_variables[1], self._start_variables[1])]
        return code

# class If


class BlockDescription(object):
    def __init__(self, block):
        block.copy_block_attributes(self)
        return

class Algorithm(object):
    def __init__(self, innerloop=True, multi_block=False, datatype="double", ):
        self.MultiBlock = multi_block
        self.subloop = innerloop
        return

    @property
    def block_descriptions(self):
        return self._block_descriptions
    
    @block_descriptions.setter
    def block_descriptions(self, blocks):
        if self.MultiBlock:
            for b in range(blocks.nblocks):
                self._block_descriptions += [BlockDescription(blocks.get_block(b))]
        return
    
    @property
    def defnintions_and_declarations(self):
        return self._def_decs
    
    @defnintions_and_declarations.setter
    def defnintions_and_declarations(self, blocks):
        defdecs = DefDecs()
        if self.MultiBlock:
            for block_number in range(blocks.nblocks):
                b = blocks.get_block(block_number)
                defdecs.add_components(b.constants.values())
                defdecs.add_components(b.Rational_constants.values())
                defdecs.add_components(b.block_datasets.values())
                defdecs.add_components(b.block_stencils.values())
        self._def_decs = defdecs
        return
    
    @property
    def nonsimulation_start(self):
        return self._nonsimulation_start
    
    @nonsimulation_start.setter
    def nonsimulation_start(self, kernels):
        self._nonsimulation_start = kernels
        return
    
    @property
    def nonsimulation_end(self):
        return self._nonsimulation_end
    
    @nonsimulation_end.setter
    def nonsimulation_end(self, kernels):
        self._nonsimulation_start = kernels
        return

class TraditionalAlgorithmRKMB(object):
    """ It is where the algorithm is generated, This is a seperate layer
    which gives user control to do any modifications for extra functionality that
    is to be performed like, doing some post processing for every time loop or
    sub rk loop
    """

    def __init__(self, blocks, dtype=None):
        self.block_descriptions = []
        self.ntimers = 0
        self.MultiBlock = True
        self.simulation_monitor = False
        if dtype:
            self.dtype = dtype
        else:
            self.dtype = "double"
        self.check_temporal_scheme(blocks)
        #self._spatial_kernels = dict(zip([no for no in range(blocks.nblocks]))
        self.prg = MainPrg()
        self.add_block_names(blocks)
        defdecs = self.get_definitions_declarations(blocks)
        self.defnitionsdeclarations = defdecs
        self.spatial_solution(blocks)
        # Now try the algorithm generation
        return

    def add_block_names(self, blocks):
        if self.MultiBlock:
            for b in range(blocks.nblocks):
                self.block_descriptions += [BlockDescription(blocks.get_block(b))]
        return

    def add_def_decs(self, defdecs):
        self.prg.add_components(defdecs)
        return

    def get_definitions_declarations(self, blocks):
        defdecs = DefDecs()
        if self.MultiBlock:
            for block_number in range(blocks.nblocks):
                b = blocks.get_block(block_number)
                defdecs.add_components(b.constants.values())
                defdecs.add_components(b.Rational_constants.values())
                defdecs.add_components(b.block_datasets.values())
                defdecs.add_components(b.block_stencils.values())
        return defdecs

    def comapre_no_sims(self, s1, s2):
        return cmp(s1.order, s2.order)
    
    #@property
    #def spatial_kernels(self):
        #return self._spatial_kernels
    
    #@spatial_kernels.setter
    #def spatial_kernels(self, block_no, kernels):
        #self._spatial_kernels[block_no] += kernels
        #return
    

    def spatial_solution(self, blocks):
        """ Add the spatial kernels to the temporal solution i.e temporalscheme.solution
        """
        print("Writing algorithm")
        fname = 'algorithm.tex'
        latex = LatexWriter()
        latex.open(fname, "Algorithm for the equations")
        if self.MultiBlock:
            bc_kernels = []
            inner_temporal_advance_kernels = []
            temporal_start = []
            temporal_end = []
            spatial_kernels = []
            before_time = []
            after_time = []
            in_time = []
            non_simulation_eqs = []
            metrics = []
            tloop_blocks = []
            inner_loop_blocks = []
            for block_number in range(blocks.nblocks):
                b = blocks.get_block(block_number)
                for scheme in b.get_temporal_schemes:
                    #print scheme
                    inner_loop_blocks += [scheme.stage]
                    tloop_blocks += [scheme.temporal_iteration]
                    for key, value in scheme.solution.iteritems():
                        #print key
                        if isinstance(key, SimulationEquations):
                            # Solution advancement kernels
                            temporal_start += scheme.solution[key].start_kernels
                            temporal_end += scheme.solution[key].end_kernels
                            inner_temporal_advance_kernels += scheme.solution[key].kernels
                            bc_kernels += key.boundary_kernels
                            spatial_kernels += key.all_spatial_kernels(b)
                        elif isinstance(key, MetricsEquation):
                            metrics += [key]
                        elif isinstance(key, NonSimulationEquations):  # Add all other types of equations
                            non_simulation_eqs += [key]
                        else:
                            if not isinstance(key, ConstituentRelations):
                                print("NOT classified", type(key))
                                raise ValueError("Equations class can not be classified: %s" % key)
            # Place the non-simulation equation kernels in the appropriate position in the simulation code
            for key in sorted(non_simulation_eqs, key=lambda x: x.order):
                for place in key.algorithm_place:
                    if isinstance(place, BeforeSimulationStarts):
                        before_time += key.Kernels
                    elif isinstance(place, AfterSimulationEnds):
                        after_time += key.Kernels
                    else:
                        if place.frequency:
                            raise NotImplementedError("In Non-simulation equations")
                        else:
                            in_time += key.Kernels
            # Process the metircs we will control here it self later we will move this to multi block
            # The first derivatives this includes bc application
            for m in metrics:
                before_time += m.fd_kernels
            # Now the second derivatives
            for m in metrics:
                before_time += m.sd_kernels
            #tloop_blocks = set(tloop_blocks)
            #inner_loop_blocks = set(inner_loop_blocks)
            if len(set(tloop_blocks)) != 1 or len(set(inner_loop_blocks)) != 1:
                raise ValueError("")
            # RK time stepping loop
            innerloop = DoLoop(inner_loop_blocks[0])
            innerloop.add_components(spatial_kernels + inner_temporal_advance_kernels + bc_kernels)

            print("Found %d kernels." % len(spatial_kernels + inner_temporal_advance_kernels + bc_kernels))
            ## Add BC kernels to temporal start
            temporal_start = bc_kernels + temporal_start
            #temporal_iteration = sc.temporal_iteration
            tloop = DoLoop(tloop_blocks[0])
            self.get_io_mb(blocks, before_time, in_time, after_time, tloop_blocks[0])
            tloop.add_components(temporal_start)
            tloop.add_components(innerloop)
            tloop.add_components(in_time)
            tloop.add_components(temporal_end)
            # Create a time loop and add the components to it
            timed_tloop = self.add_timers(tloop)
            self.prg.add_components(before_time)
            self.prg.add_components(timed_tloop)
            self.prg.add_components(after_time)
        latex.close()
        return
    
    def get_io_mb(self, blocks, before_time, in_time, after_time, temporal_iteration):
        # TODO
        for block_number in range(blocks.nblocks):
            b = blocks.get_block(block_number)
            for io in b.InputOutput:
                io.block_number = block_number # MBCHANGE
                # Check all of the datasets added to the IO classes have been defined elsewhere in the simulation
                io.check_datasets(b)
                for place in io.algorithm_place:
                    if isinstance(place, BeforeSimulationStarts):
                        io_copy = copy.deepcopy(io)
                        io_copy.dynamic_fname = False
                        before_time += [io_copy]
                    elif isinstance(place, AfterSimulationEnds):
                        io_copy = copy.deepcopy(io)
                        io_copy.dynamic_fname = False
                        after_time += [io_copy]
                    elif isinstance(place, InTheSimulation):
                        t = (Equality((temporal_iteration + 1) % place.frequency, 0))
                        cond = Condition(t)
                        io_copy = copy.deepcopy(io)
                        io_copy.dynamic_fname = True
                        io_copy.control_parameter = temporal_iteration + 1
                        cond.add_components(io_copy)
                        in_time += [cond]
                    else:
                        raise NotImplementedError("In Nonsimulation equations")
        return 

    def add_timers(self, components):
        timer = Timers(self.ntimers)
        self.ntimers += 1
        timer.add_components(components)
        return timer

    def check_temporal_scheme(self, blocks):
        """
        If Multi-block this checks the temporal scheme is the same for all the blocks
        """
        if self.MultiBlock:
            sc = set()
            for b in range(blocks.nblocks):
                for s in blocks.get_block(b).get_temporal_schemes:
                    sc.add(s.__class__.__name__)
            if len(sc) != 1:
                raise ValueError("More than one temporal scheme")
        else:
            if len(blocks[0].get_temporal_schemes) > 1:
                raise ValueError("More than one temporal scheme for a block")
        return

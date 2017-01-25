class Loop(object):
    pass

class StartLoop(Loop):
    """ Dummy place holder for various loop start"""
    pass
class EndLoop(Loop):
    """ Dummy place holder for various loop ends"""
    pass
class Doloopst(Loop):
    def __init__(self, index, name, ranges):
        self.doloopname = "%s"%name
        self.components = []
        self.loop = Idx(str(name), tuple(ranges))
        return
    def __str__(self):
        return "%s_%s"%(self.__class__.__name__, self.doloopname)
    def add_components(self, components):
        self.components += components
        return
    def get_components(self):
        return
    def opsc(self):
        return "for(int %s=%s, %s<%s, %s++)\n{"%(self.loop, str(self.loop.lower), self.loop, str(self.loop.upper), self.loop)


class Doloopend(Loop):
    def __init__(self, index):
        self.doloopname = "doloop%d"%index
        self.components = []
        return
    def __str__(self):
        return "%s_%s"%(self.__class__.__name__, self.doloopname)
    def add_components(self, components):
        self.components = []
        return
    def get_components(self):
        return
    def opsc(self):
        return "}"

class ProgramStart(StartLoop):
    """ Dummy place holder for main program, as of now it is not used but later we will move the
    logic of main program definitions to here from OPSC"""
    def __init__(self):
        return
    def name(self):
        self.name = "%sNode"%(type(self).__name__)
        return
    def components(self):

        return

class ProgramEnd(EndLoop):
    """ Dummy place holder for main program, as of now it is not used but later we will move the
    logic of main program definitions to here from OPSC"""
    def __init__(self):
        return
    def name(self):
        self.name = "%sNode"%(type(self).__name__)
        return

class PreProcess():
    def __init__(self):
        self.computations = []
        return
    def set_computations(self, computations):
        self.computations += flatten([computations])
        return
    def name(self):
        self.name = "%sNode"%(type(self).__name__)
        return

class MainPrg():
    def __init__(self, *components):
        self.components = []
        return
    def __str__(self):
        return "%s"%(self.__class__.__name__)
    def add_components(self, components):
        """ """
        self.components += components
        return
    def get_components(self):
        return

    def insert_component(self, at, components):
        location = self.locate_component(at)
        self.components[location].add_components(components)
        return

    def locate_component(self, comp):
        """Helper function to locate the component in the main program
        May be alter need to add traversing the loop to locate the component of components"""
        for no,c in enumerate(self.components):
            if str(c) == str(comp):
                return no
from .block import SimulationBlock as SB
from .kernel import Kernel
from .latex import *
class TraditionalAlgorithm(MainPrg):
    def __init__(self, blocks):
        if isinstance(blocks, SB):
            blocks = [blocks]
        self.generate_time_loop(blocks)
        return

    def generate_time_loop(self, blocks):
        """ Add the spatial kernels to the temporal solution now the scheme.solution
        """
        fname = './algorithm.tex'
        latex = LatexWriter()
        latex.open('./algorithm.tex')
        metadata = {"title": "Algorithm for the equations", "author": "Jammy", "institution": ""}
        latex.write_header(metadata)
        for b in blocks:
            for scheme in b.get_temporal_schemes:
                for key, value in scheme.solution.iteritems():
                    for l in b.list_of_equation_classes:
                        if l.order >=0 and l.order <100: #Checks if the equation classes are part of the time loop
                            scheme.solution[key].kernels = l.all_spatial_kernels + scheme.solution[key].kernels
                        else:
                            NotImplementedError("")
                #print [s.computation_name for s in scheme.solution[key].kernels]
                for s in scheme.solution[key].kernels:
                    s.write_latex(latex)
        latex.write_footer()
        latex.close()
        return


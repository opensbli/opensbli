from sympy import flatten, Max
from .latex import *
class Kernel(object):

    """ A computational kernel which will be executed over all the grid points and in parallel. """
    mulfactor = {0:-1, 1:1}
    def __init__(self, block, computation_name = None):
        """ Set up the computational kernel"""
        self.block_number = block.blocknumber
        self.computation_name = computation_name
        #self.kernel_number = block.kernel_counter
        #block.increase_kernel_counter
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
        from sympy.tensor import IndexedBase
        from .opensbliobjects import DataSet
        requires = []
        for eq in self.equations:
            requires += list(eq.rhs.atoms(DataSet))
        #print requires
        return requires

    def write_latex(self, latex):
        latex.write_string('The kernel is %s'%self.computation_name)
        #latex.write_string('The halos are ')
        for index, eq in enumerate(self.equations):
            if isinstance(eq, Equality):
                latex.write_expression(eq)
        return

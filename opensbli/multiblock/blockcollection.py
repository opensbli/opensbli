
from opensbli.core.block import SimulationBlock
from opensbli.core.boundary_conditions.multi_block import MultiBlockBoundary
import copy
from sympy import flatten

# MBCHANGE
@property
def get_interface_bc(self):
    interface_bcs = []
    for no, b in enumerate(self.boundary_types):
        for side in [0,1]:
            bc = self.boundary_types[no][side]
            if isinstance(bc, MultiBlockBoundary):
                interface_bcs += [bc]
    return interface_bcs

def apply_interface_bc(self, arrays, multiblock_description):
    kernels = []
    interface_bcs = self.get_interface_bc
    for bc in interface_bcs:
        k = bc.apply_interface(arrays, self, multiblock_description)
        if isinstance(k, list):
            kernels += k
        else:
            kernels += [k]
    return kernels
# MBCHANGE
SimulationBlock.get_interface_bc = get_interface_bc
SimulationBlock.apply_interface_bc = apply_interface_bc
# MBCHANGE

class MultiBlock():
    def __init__(self, ndim, nblocks):
        self.nblocks = nblocks
        self.blocks = [0 for i in range(nblocks)]
        for i in range(nblocks):
            self.blocks[i] = SimulationBlock(ndim, block_number=i)
        return

    def get_block(self, number):
        return self.blocks[number]

    def set_discretisation_schemes(self, schemes):
        for i in range(self.nblocks):
            schemes1 = copy.deepcopy(schemes)
            self.blocks[i].set_discretisation_schemes(schemes1)
        return

    def set_equations(self, list_of_equations):
        for b in self.blocks:
            b.set_equations([copy.deepcopy(e) for e in list_of_equations])
        return

    def set_filters(self, list_of_filters):
        for i, b in enumerate(self.blocks):
            filt = flatten(list_of_filters[i])
            b.set_equations(filt)
        return
    
    def set_block_boundaries(self, bclist):    
        if len(bclist.keys()) != self.nblocks:
            raise ValueError("")
        for b in self.blocks:
            b.set_block_boundaries(bclist[b.blocknumber])
        return
    
    def discretise(self):
        for b in self.blocks:
            b.discretise()
        # After discretisation apply the interface boundary conditions as this requires halos required to be populated
        # for the other blocks
        self.apply_interface_bc()
        return
    def apply_interface_bc(self):
        for block in self.blocks:
            for eq in block.list_of_equation_classes:
                eq.apply_interface_bc(block, self)
        return
            
    def setio(self, list_of_ios):
        for b in self.blocks:
            copied_io = [copy.deepcopy(io) for io in list_of_ios]
            b.setio(copied_io)
        return
# MBCHANGE
from sympy import flatten
side_names = {0: 'left', 1: 'right'}


class Exchange(object):
    pass


class ExchangeSelf(Exchange):
    """ Defines data exchange on the same block. """

    def __init__(self, block, direction, side):
        # Range of evaluation (i.e. the grid points, including the halo points, over which the computation should be performed).
        self.computation_name = "exchange"
        self.block_number = block.blocknumber
        self.block_name = block.blockname
        self.direction = direction
        self.side = side_names[side]
        self.flip = [False]  # MBCHANGE
        return

    @property
    def name(self):
        return "%s%d_block%d" % (self.computation_name, self.number, self.block_number)  # MBCHANGE

    def set_arrays(self, arrays):
        self.transfer_arrays = flatten(arrays)
        self.from_arrays = flatten(arrays)
        self.to_arrays = flatten(arrays)
        return

    def set_transfer_from(self, transfer):
        self.transfer_from = transfer
        return

    def set_transfer_to(self, transfer):
        self.transfer_to = transfer
        return

    def set_transfer_size(self, size):
        self.transfer_size = size
        return

    @property
    def algorithm_node_name(self):
        name = "Boundary_exchange_block_%d_direction%d_side%s" % (self.block_number, self.direction, self.side)
        return name

    def write_latex(self, latex):
        latex.write_string("This is an exchange self kernel on variables %s\\\\" % ', '.join([str(a) for a in self.transfer_arrays]))
        return

    @property
    def opsc_code(self):
        """The string for calling the boundary condition in OPSC is updated while creating
        the code for exchanging data."""
        return [self.call_name]

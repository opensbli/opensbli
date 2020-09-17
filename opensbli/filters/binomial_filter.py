""" David J. Lusher: Binomial filter using the coefficients of (a-b)^n / 2^n for filter of order n.
    Original three-point 2nd order version by Alex Gillespie."""

from opensbli import *
from sympy import symbols, exp, pprint, Piecewise, binomial
from opensbli.core.opensbliobjects import DataObject, ConstantObject, GroupedPiecewise
from opensbli.equation_types.opensbliequations import OpenSBLIEquation
from opensbli.postprocess.post_process_eq import *
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.code_generation.algorithm.common import *
from opensbli.utilities.user_defined_kernels import UserDefinedEquations
from opensbli.core.block import SimulationBlock
from opensbli.multiblock.blockcollection import MultiBlock

class BinomialFilter(object):
    def __init__(self, block, order, grid_condition=None, sigma=0.05):
        self.filter_no = block.blocknumber
        if (order % 2) != 0:
            raise ValueError("The filter is only defined for even orders n.")
        elif (order > 10):
            raise ValueError("Increase the number of halo points in scheme.py for high order filters.")
        else:
            self.order = order
        # Spatial dependence of the filter
        self.grid_condition = grid_condition
        # Width and weightings of the filter
        self.generate_weights()
        sigma_symbol = ConstantObject('sigma_filt')
        sigma_symbol.value = sigma
        self.sigma = sigma_symbol
        self.equation_classes = []
        # Create the filter equations
        self.create_filter(block)
        return

    def generate_weights(self):
        """ Creates the binomial coefficients for a filter of order N."""
        N = self.order
        self.weights = [binomial(N, i)/2.0**N for i in range(N + 1)]
        self.locations = [i for i in range(-int(N/2.0), int(N/2.0)+1)]
        print("Using a binomial filter of order %d for block %d." % (N, self.filter_no))
        return

    def create_stencil(self, block, q, direction):
        """ Indexes the datasets based on the width of the filter stencil."""
        output = []
        for dset in q:
            dset = block.location_dataset(dset)
            dset_locations = []
            for i, location in enumerate(self.locations):
                dset_locations.append(increment_dataset(dset, direction, location))
            output += [dset_locations]
        return output

    def filtered_equations(self, q_f, q_stencil):
        """ Creates the weighted sum for the binomial filter."""
        output = []
        for u_f, u_stencil in zip(q_f, q_stencil):
            rhs = [self.weights[i]*u_stencil[i] for i, location in enumerate(self.locations)]
            output += [OpenSBLIEquation(u_f, sum(rhs))]
        return output

    def conditional_expression(self, filter_equations):
        """ Only applies the filter in certain regions of the domain. Based on the 
        grid_condition input to the class, which should be a boolean expression built from the coordinate arrays."""
        output_equations = []
        for i, eqn in enumerate(filter_equations):
            condition = ExprCondPair(eqn.rhs, self.grid_condition)
            output_equations += [OpenSBLIEquation(eqn.lhs, Piecewise(condition, (eqn.lhs, True)))]
        return output_equations

    def create_equations(self, block):
        ndim = block.ndim
        # Conservative variables
        if ndim == 2:
            q = ['rho', 'rhou0', 'rhou1', 'rhoE']
        elif ndim == 3:
            q = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']

        # Create the three point stencils
        q_xstencil = self.create_stencil(block, q, 0)
        q_ystencil = self.create_stencil(block, q, 1)
        q_fx = [GridVariable(u + "_xfiltered") for u in q]
        q_fy = [GridVariable(u + "_yfiltered") for u in q]
        if ndim == 3:
            q_zstencil = self.create_stencil(block, q, 2)
            q_fz = [GridVariable(u + "_zfiltered") for u in q]
        
        # Create the filter equations
        output_equations = self.filtered_equations(q_fx, q_xstencil)
        output_equations += self.filtered_equations(q_fy, q_ystencil)
        if ndim == 3:
            output_equations += self.filtered_equations(q_fz, q_zstencil)
        # Average the filter
        q_f = [GridVariable(u + "_filtered") for u in q]
        if ndim == 2:
            for u_f, u_fx, u_fy in zip(q_f, q_fx, q_fy):
                output_equations += [OpenSBLIEquation(u_f, (u_fx + u_fy)/2.0)]
        elif ndim == 3:
            for u_f, u_fx, u_fy, u_fz in zip(q_f, q_fx, q_fy, q_fz):
                output_equations += [OpenSBLIEquation(u_f, (u_fx + u_fy + u_fz)/3.0)]
        # Blend the filter
        blended_equations = []
        for u, u_f in zip(q, q_f):
            u = block.location_dataset(u)
            blended_equations += [OpenSBLIEquation(u, (1-self.sigma)*u + self.sigma*u_f)]
        # Apply the spatial condition
        if self.grid_condition is not None:
            output_equations += self.conditional_expression(blended_equations)
        else:
            output_equations += blended_equations
        return output_equations

    def create_filter(self, block):
        # Create a kernel at the end of the time loop, every iteration (no frequency)
        filter_class = UserDefinedEquations()
        filter_class.algorithm_place = InTheSimulation(frequency=False)
        filter_class.computation_name = 'Binomial filter'
        # Place the filter at the very end
        filter_class.order = 10000
        filter_class.add_equations(self.create_equations(block))
        self.equation_classes.append(filter_class)
        return

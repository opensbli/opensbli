
#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (c) see License file

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

"""@brief
   @authors David J Lusher
   @contributors
   @details
"""

from sympy import IndexedBase, Symbol, Rational, Abs, flatten, Max, S, ceiling
from opensbli.core.opensblifunctions import TenoDerivative
from opensbli.core.opensbliobjects import ConstantObject
from opensbli.core.grid import GridVariable
from opensbli.schemes.spatial.scheme import Scheme
from opensbli.schemes.spatial.weno import LLFCharacteristic, ShockCapturing, RFCharacteristic
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.equation_types.opensbliequations import OpenSBLIEq, SimulationEquations


class TenoHalos(object):
    """ Object for TENO halos.

    :arg int order: Order of the TENO scheme.
    :arg bool reconstruction: True if halos for a reconstruction. """

    def __init__(self, order, reconstruction=None):
        if not reconstruction:
            if order == 8:
                k = 4
            else:
                k = 3  # Halos correct for TENO5, TENO6
            self.halos = [-k, k+1]
        else:
            self.halos = [-1, 1]
        return

    def get_halos(self, side):
        return self.halos[side]


class TenoStencil(object):
    """ Stencil object used for the TENO reconstruction.

    :arg int side: Side of the TENO reconstruction, either -1 (left) or +1 (right).
    :arg int width: Width of the TENO candidate stencil.
    :arg int stencil_number: Index of the stencil. """

    def __init__(self, side, width, stencil_number):
        self.width = width
        self.side = side
        self.stencil_number = stencil_number
        self.fn_points = None
        self.eno_coeffs = None
        self.opt_coeff = None
        self.smoothness_indicator = None
        return


class ConfigureTeno(object):
    """ Object containing the parameters needed by the TENO reconstruction for a given order and side.

    :arg int order: Order of the TENO scheme.
    :arg int side: Side of the TENO reconstruction, either -1 (left) or +1 (right).
    :arg bool optimized: Optimized or regular coefficients. """

    def __init__(self, order, side, optimized=False):
        self.order = order
        self.side = side
        self.optimized = optimized
        self.n_stencils = order - 2
        # Create stencils
        self.stencils = self.generate_stencils()
        self.fn_points = self.generate_func_points()
        self.eno_coeffs = self.generate_eno_coefficients()
        self.opt_coeffs = self.generate_optimal_coefficients()
        self.fn_dictionary, self.smoothness_indicators, self.smoothness_symbols = self.generate_smoothness_indicators()
        self.unique_fn_points = sorted(set(flatten(self.fn_points)))
        # Add the function points and coefficients to the stencil objects
        self.update_stencils()
        return

    def generate_stencils(self):
        """ Create stencils required for TENO5, TENO6 and TENO8.

        :return: stencils: List of TENO stencil objects.
        :rtype: list"""
        stencils = []
        side = self.side
        if self.order == 5:
            widths = [3, 3, 3]
        elif self.order == 6:
            widths = [3, 3, 3, 4]
        elif self.order == 8:
            widths = [3, 3, 3, 4, 4, 5]
        for i, width in enumerate(widths):
            stencils += [TenoStencil(side, width, i)]
        return stencils

    def generate_func_points(self):
        """ Generates the relative function points required for the TENO stencils.

        :returns: fn_points: List of stencil point lists."""
        fn_points = []
        side, order = self.side, self.order
        if side == 1:
            fn_points = [[-1, 0, 1], [0, 1, 2], [-2, -1, 0], [0, 1, 2, 3], [-3, -2, -1, 0], [0, 1, 2, 3, 4]]
        elif side == -1:
            fn_points = [[0, 1, 2], [-1, 0, 1], [1, 2, 3], [-2, -1, 0, 1], [1, 2, 3, 4], [-3, -2, -1, 0, 1]]
        return [fn_points[i] for i in range(order-2)]

    def generate_eno_coefficients(self):
        """ Generates the ENO coefficients required for the TENO stencils.

        :returns: coeffs: List of ENO coefficients for each stencil."""
        side = self.side
        if side == 1:
            coeffs = [[Rational(-1, 6), Rational(5, 6), Rational(2, 6)], [Rational(2, 6), Rational(5, 6), Rational(-1, 6)],
                      [Rational(2, 6), Rational(-7, 6), Rational(11, 6)]]
            if self.order in [6, 8]:
                coeffs += [[Rational(3, 12), Rational(13, 12), Rational(-5, 12), Rational(1, 12)]]
        elif side == -1:
            coeffs = [[Rational(2, 6), Rational(5, 6), Rational(-1, 6)], [Rational(-1, 6), Rational(5, 6), Rational(2, 6)],
                      [Rational(11, 6), Rational(-7, 6), Rational(1, 3)]]
            if self.order in [6, 8]:
                coeffs += [[Rational(1, 12), Rational(-5, 12), Rational(13, 12), Rational(3, 12)]]
        return coeffs

    def generate_optimal_coefficients(self):
        """ Generates the optimal linear coefficients required for the TENO stencils.

        :returns: opt_coeffs: List of optimal coefficients for each stencil."""
        # Optimal weights are not reversed for TENO when switching between upwind/downwind biasing
        order = self.order
        if self.optimized:
            if order == 5:
                opt_coeffs = [Rational(11, 20), Rational(4, 10), Rational(1, 20)]
            elif order == 6:
                opt_coeffs = [Rational(231, 500), Rational(3, 10), Rational(27, 500), Rational(23, 125)]
        elif not self.optimized:
            if order == 5:
                opt_coeffs = [Rational(6, 10), Rational(3, 10), Rational(1, 10)]
            elif order == 6:
                opt_coeffs = [Rational(9, 20), Rational(6, 20), Rational(1, 20), Rational(4, 20)]
        return opt_coeffs

    def generate_smoothness_indicators(self):
        """ Generates the smoothness indicators required for the TENO stencils.

        :returns: fns_dictionary: Key: Integer grid location, Value: placeholder function 'f' at the grid location.
        :returns: smoothness_indicators: List of smoothness indicator expressions.
        :returns: smoothness_symbols: List of placeholder symbols for the smoothness indicators."""
        points = sorted(list(set(flatten([i for i in self.fn_points]))))
        f = IndexedBase('f')
        symbolic_functions = []
        fns = {}
        for p in points:
            symbolic_functions.append(f[p])
            fns[p] = symbolic_functions[-1]
        smoothness_indicators = []
        if self.side == 1:
            # Stencils [0, 1, 2] for TENO5
            # Factored versions: 0, 1, 2, 3
            smoothness_indicators += [Rational(1, 4)*(fns[-1] - fns[1])**2 + Rational(13, 12)*(fns[-1]-2*fns[0]+fns[1])**2]
            smoothness_indicators += [Rational(1, 4)*(3*fns[0]-4*fns[1]+fns[2])**2 + Rational(13, 12)*(fns[0]-2*fns[1]+fns[2])**2]
            smoothness_indicators += [Rational(1, 4)*(fns[-2]-4*fns[-1]+3*fns[0])**2 + Rational(13, 12)*(fns[-2]-2*fns[-1]+fns[0])**2]
            if self.order in [6, 8]:
                smoothness_indicators += [Rational(1, 36)*(-11*fns[0]+18*fns[1]-9*fns[2]+2*fns[3])**2 + Rational(13, 12)*(2*fns[0]-5*fns[1]+4*fns[2]-fns[3])**2 +
                                          Rational(781, 720)*(-fns[0]+3*fns[1]-3*fns[2]+fns[3])]
        elif self.side == -1:
            # Stencils [0, 1, 2] for TENO5
            # Factored versions: 0, 1, 2, 3 [0, 1, 2], [-1, 0, 1], [1, 2, 3], [-2, -1, 0, 1]
            smoothness_indicators += [Rational(1, 4)*(fns[0]-fns[2])**2 + Rational(13, 12)*(fns[0]-2*fns[1]+fns[2])**2]
            smoothness_indicators += [Rational(1, 4)*(fns[-1]-4*fns[0]+3*fns[1])**2 + Rational(13, 12)*(fns[-1]-2*fns[0]+fns[1])**2]
            smoothness_indicators += [Rational(1, 4)*(3*fns[1]-4*fns[2]+fns[3])**2 + Rational(13, 12)*(fns[1]-2*fns[2]+fns[3])**2]
            if self.order in [6, 8]:
                smoothness_indicators += [Rational(1, 36)*(-2*fns[-2]+9*fns[-1]-18*fns[0]+11*fns[1])**2 +
                                          Rational(13, 12)*(-1*fns[-2]+4*fns[-1]-5*fns[0]+2*fns[1])**2 +
                                          Rational(781, 720)*(-1*fns[-2]+3*fns[-1]-3*fns[0]+fns[1])**2]
        fns_dictionary = fns
        smoothness_symbols = [Symbol('beta_%d' % r) for r in range(len(smoothness_indicators))]
        return fns_dictionary, smoothness_indicators, smoothness_symbols

    def update_stencils(self):
        """Function to update the TenoStencil objects with their respective properties.

        :returns: None """
        for stencil in self.stencils:
            no = stencil.stencil_number
            stencil.fn_points = self.fn_points[no]
            stencil.eno_coeffs = self.eno_coeffs[no]
            stencil.opt_coeff = self.opt_coeffs[no]
            stencil.smoothness_indicator = self.smoothness_indicators[no]
        return


class Teno5(object):
    """ Base class for 5th order TENO scheme."""

    def __init__(self):
        return

    def generate_alphas(self, RV, TC):
        """ Create the alpha terms for the non-linear TENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object TC: Configuration settings for a reconstruction of either left or right."""
        # Scale separation parameters
        C, q = S.One, 6
        # Global reference smoothness indicator tau_5
        tau_5 = Abs(RV.smoothness_symbols[0] - RV.smoothness_symbols[-1])
        for r in range(TC.n_stencils):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append((C + (tau_5/(self.eps + RV.smoothness_symbols[r])))**q)
        return


class Teno6(object):
    """ Base class for 6th order TENO scheme."""

    def __init__(self):
        return

    def generate_alphas(self, RV, TC):
        """ Create the alpha terms for the non-linear TENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object TC: Configuration settings for a reconstruction of either left or right."""
        # Scale separation parameters
        C, q = S.One, 6
        # Global reference smoothness indicator tau_6
        tau_6 = RV.smoothness_symbols[3] - Rational(1, 6)*(RV.smoothness_symbols[0] + RV.smoothness_symbols[2] - 4*RV.smoothness_symbols[1])
        for r in range(TC.n_stencils):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append((C + (tau_6/(self.eps + RV.smoothness_symbols[r])))**q)
        return


class TenoReconstructionVariable(object):
    """ Reconstruction variable object to hold the quantities required for TENO.

    :arg str name: Name of the reconstruction, either left or right."""

    def __init__(self, name):
        self.name = name
        self.smoothness_indicators = []
        self.smoothness_symbols = []
        self.alpha_evaluated = []
        self.alpha_symbols = []
        self.inv_alpha_sum_symbols = []
        self.inv_alpha_sum_evaluated = []
        self.inv_omega_sum_symbols = []
        self.inv_omega_sum_evaluated = []
        self.tau_8_symbol = []
        self.tau_8_evaluated = []
        self.omega_evaluated = []
        self.omega_symbols = []
        self.kronecker_evaluated = []
        self.kronecker_symbols = []
        self.function_stencil_dictionary = {}
        self.reconstructed_expression = None
        self.reconstructed_symbol = GridVariable('%s' % (name))
        return

    def update_quantities(self, original):
        """ Updates the quantities required by TENO in the reconstruction variable.

        :arg object original: Reconstruction object variable, either left or right reconstruction."""
        self.smoothness_symbols += [GridVariable('%s' % (s)) for s in original.smoothness_symbols]
        self.alpha_symbols += [GridVariable('%s' % (s)) for s in original.alpha_symbols]
        self.inv_alpha_sum_symbols += [GridVariable('%s' % (s)) for s in original.inv_alpha_sum_symbols]
        self.inv_omega_sum_symbols += [GridVariable('%s' % (s)) for s in original.inv_omega_sum_symbols]
        self.omega_symbols += [GridVariable('%s' % (s)) for s in original.omega_symbols]
        self.kronecker_symbols += [GridVariable('%s' % (s)) for s in original.kronecker_symbols]

        if original.order == 8:
            self.tau_8_symbol = [GridVariable('%s' % (original.tau_8_symbol[0]))]
            originals = original.tau_8_symbol
            new = self.tau_8_symbol[:]
        else:
            originals = []
            new = []

        originals += original.smoothness_symbols + original.alpha_symbols + original.inv_alpha_sum_symbols + original.kronecker_symbols + original.inv_omega_sum_symbols + original.omega_symbols
        new += self.smoothness_symbols + self.alpha_symbols + self.inv_alpha_sum_symbols + self.kronecker_symbols + self.inv_omega_sum_symbols + self.omega_symbols
        subs_dict = dict(zip(originals, new))

        for key, value in original.function_stencil_dictionary.iteritems():
            subs_dict[value] = self.function_stencil_dictionary[key]

        self.smoothness_indicators = [s.subs(subs_dict) for s in original.smoothness_indicators]
        # Only update tau_8 for right reconstructions
        if original.order == 8:
            if isinstance(self, LeftTenoReconstructionVariable):
                self.tau_8_evaluated = self.tau_8_symbol
            else:
                self.tau_8_evaluated = [s.subs(subs_dict) for s in original.tau_8_evaluated]

        self.alpha_evaluated = [s.subs(subs_dict) for s in original.alpha_evaluated]
        self.inv_alpha_sum_evaluated = [s.subs(subs_dict) for s in original.inv_alpha_sum_evaluated]
        self.inv_omega_sum_evaluated = [s.subs(subs_dict) for s in original.inv_omega_sum_evaluated]
        self.omega_evaluated = [s.subs(subs_dict) for s in original.omega_evaluated]
        self.kronecker_evaluated = [s.subs(subs_dict) for s in original.kronecker_evaluated]
        self.reconstructed_expression = original.reconstructed_expression.subs(subs_dict)
        return

    def evaluate_quantities(self):
        """ Adds the evaluations of the TENO quantities to the computational kernel.

        :arg object kernel: OpenSBLI Kernel for the TENO computation."""
        all_symbols = self.smoothness_symbols + self.tau_8_symbol + self.alpha_symbols + self.inv_alpha_sum_symbols + self.kronecker_symbols + self.inv_omega_sum_symbols + self.omega_symbols
        all_evaluations = self.smoothness_indicators + self.tau_8_evaluated + self.alpha_evaluated + self.inv_alpha_sum_evaluated + self.kronecker_evaluated + self.inv_omega_sum_evaluated + self.omega_evaluated

        final_equations = []
        for no, value in enumerate(all_symbols):
            final_equations += [OpenSBLIEq(value, all_evaluations[no], evaluate=False)]
        self.final_equations = final_equations
        rv = self.reconstructed_symbol
        if self.settings.has_key("combine_reconstructions") and self.settings["combine_reconstructions"]:
            self.final_equations += [OpenSBLIEq(rv, rv + self.reconstructed_expression)]
        else:
            self.final_equations += [OpenSBLIEq(rv, self.reconstructed_expression)]
        return


class LeftTenoReconstructionVariable(TenoReconstructionVariable):
    """ Reconstruction object for the left TENO reconstruction.

        :arg str name: 'left' """

    def __init__(self, name):
        TenoReconstructionVariable.__init__(self, name)
        return


class RightTenoReconstructionVariable(TenoReconstructionVariable):
    """ Reconstruction object for the right TENO reconstruction.

        :arg str name: 'right' """

    def __init__(self, name):
        TenoReconstructionVariable.__init__(self, name)
        return


class Teno(Scheme, ShockCapturing):
    """ The main TENO class, part of the ShockCapturing family.

        :arg int order: Numerical order of the TENO scheme (3,5,...)."""

    def __init__(self, order, **kwargs):
        Scheme.__init__(self, "TenoDerivative", order)
        self.schemetype = "Spatial"
        self.order = order
        if order == 5:
            WT = Teno5()
        elif order == 6:
            WT = Teno6()
        else:
            raise NotImplementedError("Only 5th, 6th and 8th order TENO are currently implemented.")
        # Epsilon to avoid division by zero in non-linear weights
        WT.eps = ConstantObject('eps')
        # Parameter to control the spectral properties of the TENO scheme
        self.CT = ConstantObject('TENO_CT')
        CTD.add_constant(self.CT)
        CTD.add_constant(WT.eps)
        # Use optimized schemes?
        optimized = True
        self.halotype = TenoHalos(self.order)
        self.required_constituent_relations_symbols = {}
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        self.reconstruction_classes = [RightTenoReconstructionVariable('right'), LeftTenoReconstructionVariable('left')]
        # Populate the quantities required by TENO for the left and right reconstruction variable.
        for no, side in enumerate([1, -1]):
            TC = ConfigureTeno(order, side, optimized)
            RV = self.reconstruction_classes[no]
            RV.order = order
            RV.func_points, RV.n_stencils = sorted(set(flatten(TC.fn_points))), TC.n_stencils
            RV.stencil_points, RV.function_stencil_dictionary = TC.unique_fn_points, TC.fn_dictionary
            RV.smoothness_symbols, RV.smoothness_evaluated = TC.smoothness_symbols, TC.smoothness_indicators
            if (self.order == 8):
                WT.global_smoothness_indicator(RV)
            # Only the generation of alphas changes for different orders of TENO
            WT.generate_alphas(RV, TC)
            self.create_cutoff_equations(RV, TC)
            self.generate_omegas(RV, TC)
            self.generate_reconstruction(RV, TC)
            self.reconstruction_classes[no] = RV
        return

    def reconstruction_halotype(self, order, reconstruction=True):
        return TenoHalos(order, reconstruction)

    def group_by_direction(self, eqs):
        """ Groups the input equations by the direction (x0, x1, ...) they depend upon.

        :arg list eqs: List of equations to group by direction.
        :returns: grouped: Dictionary of {direction: equations} key, value pairs for equations grouped by direction."""
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(TenoDerivative))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped

    def generate_reconstruction(self, RV, TC):
        """ Create the final TENO stencil by summing the stencil points, ENO coefficients and TENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object TC: Configuration settings for a reconstruction of either left or right."""
        reconstruction = 0
        fns = []
        for stencil in TC.stencils:
            RV.smoothness_indicators.append(stencil.smoothness_indicator)
            fns = [RV.function_stencil_dictionary[i] for i in stencil.fn_points]
            eno_interpolation = sum([point*coefficient for (point, coefficient) in zip(fns, stencil.eno_coeffs)])
            reconstruction += RV.omega_symbols[stencil.stencil_number]*eno_interpolation
        RV.reconstructed_expression = reconstruction
        return

    def generate_omegas(self, RV, TC):
        """ Create the omega terms for the non-linear TENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object TC: Configuration settings for a reconstruction of either left or right."""
        inv_omega_sum_symbols = [Symbol('inv_omega_sum')]
        inv_omega_sum_evaluated = [S.One/sum([TC.opt_coeffs[i]*RV.kronecker_symbols[i] for i in range(RV.n_stencils)])]
        RV.omega_symbols = [Symbol('omega_%d' % r) for r in range(RV.n_stencils)]
        RV.omega_evaluated = [TC.opt_coeffs[r]*RV.kronecker_symbols[r]*inv_omega_sum_symbols[0] for r in range(RV.n_stencils)]
        RV.inv_omega_sum_symbols = inv_omega_sum_symbols
        RV.inv_omega_sum_evaluated = inv_omega_sum_evaluated
        return

    def create_cutoff_equations(self, RV, TC):
        """ Generates the discrete cut-off functions that determine whether a certain TENO stencil is to
        contribute to the final flux reconstruction. These are evaluated with a combination of ceiling and Max functions
        to avoid having to evaluate expensive if else conditions in the generated code.

        :arg object RV: The reconstruction variable object.
        :arg object TC: Configuration settings for a reconstruction of either left or right. """
        eqns = []
        kronecker_deltas = [Symbol('delta_%d' % r) for r in range(TC.n_stencils)]
        inv_alpha_sum_symbols = [Symbol('inv_alpha_sum')]
        inv_alpha_sum_evaluated = [S.One/sum(RV.alpha_symbols)]
        for r in range(RV.n_stencils):
            eqns += [ceiling((Max(RV.alpha_symbols[r]*inv_alpha_sum_symbols[0]-self.CT, S.Zero)))]
        RV.inv_alpha_sum_symbols = inv_alpha_sum_symbols
        RV.inv_alpha_sum_evaluated = inv_alpha_sum_evaluated
        RV.kronecker_evaluated = eqns
        RV.kronecker_symbols = kronecker_deltas
        return


class LLFTeno(LLFCharacteristic, Teno):
    """ Local Lax-Friedrichs flux splitting applied to characteristic variables using a TENO scheme.

    :arg int order: Order of the WENO/TENO scheme.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging."""

    def __init__(self, order, physics=None, averaging=None):
        LLFCharacteristic.__init__(self, physics, averaging)
        print "A TENO scheme of order %s is being used for shock capturing" % str(order)
        Teno.__init__(self, order)
        return

    def discretise(self, type_of_eq, block):
        """ This is the place where the logic of vector form of equations are implemented.
        Find physical fluxes by grouping derivatives by direction --> in central, copy over
        Then the physical fluxes are transformed to characteristic space ---> a function in Characteristic
        For each f+ and f-, find f_hat of i+1/2, i-1/2, (L+R) are evaluated  ----> Function in WENO scheme, called from in here
        flux at i+1/2 evaluated -- > Function in WENO scheme
        Then WENO derivative class is instantiated with the flux at i+1/2 array --> Function in WENO scheme, called from in here
        Final derivatives are evaluated from Weno derivative class --> Using WD.discretise."""
        if isinstance(type_of_eq, SimulationEquations):
            eqs = flatten(type_of_eq.equations)
            grouped = self.group_by_direction(eqs)
            all_derivatives_evaluated_locally = []
            reconstruction_halos = self.reconstruction_halotype(self.order, reconstruction=True)
            solution_vector = flatten(type_of_eq.time_advance_arrays)

            # Instantiate eigensystems with block, but don't add metrics yet
            self.instantiate_eigensystem(block)

            for direction, derivatives in grouped.iteritems():
                # Create a work array for each component of the system
                all_derivatives_evaluated_locally += derivatives
                for no, deriv in enumerate(derivatives):
                    deriv.create_reconstruction_work_array(block)
                # Kernel for the reconstruction in this direction
                kernel = self.create_reconstruction_kernel(direction, reconstruction_halos, block)
                # Get the equations for a characteristic reconstruction
                characteristic_eqns = self.get_characteristic_equations(direction, derivatives, solution_vector, block)
                pre_process, interpolated, post_process = characteristic_eqns[0], characteristic_eqns[1], characteristic_eqns[2]
                # Add the equations to the kernel and add the kernel to SimulationEquations
                kernel.add_equation(pre_process + interpolated + post_process)
                type_of_eq.Kernels += [kernel]
            # Generate kernels for the constituent relations
            if grouped:
                constituent_relations = self.generate_constituent_relations_kernels(block)
                type_of_eq.Kernels += [self.evaluate_residuals(block, eqs, all_derivatives_evaluated_locally)]
                constituent_relations = self.check_constituent_relations(block, eqs, constituent_relations)
            return constituent_relations


class RFTeno(RFCharacteristic, Teno):
    """ Roe flux splitting applied to characteristic variables using a TENO scheme.

    :arg int order: Order of the WENO/TENO scheme.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging."""

    def __init__(self, order, physics=None, averaging=None):
        RFCharacteristic.__init__(self, physics, averaging)
        print "A TENO scheme of order %s is being used for shock capturing with Roe-flux differencing" % str(order)
        Teno.__init__(self, order)
        return

    def discretise(self, type_of_eq, block):
        """ This is the place where the logic of vector form of equations are implemented.
        Find physical fluxes by grouping derivatives by direction --> in central, copy over
        Then the physical fluxes are transformed to characteristic space ---> a function in Characteristic
        For each f+ and f-, find f_hat of i+1/2, i-1/2, (L+R) are evaluated  ----> Function in WENO scheme, called from in here
        flux at i+1/2 evaluated -- > Function in WENO scheme
        Then WENO derivative class is instantiated with the flux at i+1/2 array --> Function in WENO scheme, called from in here
        Final derivatives are evaluated from Weno derivative class --> Using WD.discretise."""
        if isinstance(type_of_eq, SimulationEquations):
            eqs = flatten(type_of_eq.equations)
            grouped = self.group_by_direction(eqs)
            all_derivatives_evaluated_locally = []
            reconstruction_halos = self.reconstruction_halotype(self.order, reconstruction=True)
            solution_vector = flatten(type_of_eq.time_advance_arrays)

            # Instantiate eigensystems with block, but don't add metrics yet
            self.instantiate_eigensystem(block)

            for direction, derivatives in grouped.iteritems():
                # Create a work array for each component of the system
                all_derivatives_evaluated_locally += derivatives
                for no, deriv in enumerate(derivatives):
                    deriv.create_reconstruction_work_array(block)
                # Kernel for the reconstruction in this direction
                kernel = self.create_reconstruction_kernel(direction, reconstruction_halos, block)
                # Get the equations for a characteristic reconstruction
                characteristic_eqns = self.get_characteristic_equations(direction, derivatives, solution_vector, block)
                pre_process, interpolated, post_process = characteristic_eqns[0], characteristic_eqns[1], characteristic_eqns[2]
                # Add the equations to the kernel and add the kernel to SimulationEquations
                kernel.add_equation(pre_process + interpolated + post_process)
                type_of_eq.Kernels += [kernel]
            # Generate kernels for the constituent relations
            if grouped:
                constituent_relations = self.generate_constituent_relations_kernels(block)
                type_of_eq.Kernels += [self.evaluate_residuals(block, eqs, all_derivatives_evaluated_locally)]
                constituent_relations = self.check_constituent_relations(block, eqs, constituent_relations)
            return constituent_relations

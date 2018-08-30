
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
   @authors Satya Pramod Jammy, David J Lusher
   @contributors
   @details
"""

from sympy import IndexedBase, Symbol, Rational, solve, interpolating_poly, integrate, Abs, Float, flatten, S
from opensbli.core.opensblifunctions import WenoDerivative
from opensbli.equation_types.opensbliequations import SimulationEquations, OpenSBLIEq
from opensbli.core.grid import GridVariable
from .scheme import Scheme
from sympy import horner
from opensbli.schemes.spatial.shock_capturing import ShockCapturing, LLFCharacteristic, RFCharacteristic


class WenoHalos(object):
    """ Object for WENO halos.

    :arg int order: Order of the WENO scheme.
    :arg bool reconstruction: True if halos for a reconstruction. """

    def __init__(self, order, reconstruction=None):
        k = int(0.5*(order+1))
        if not reconstruction:
            self.halos = [-k, k+1]
        else:
            self.halos = [-1, 1]
        return

    def get_halos(self, side):
        return self.halos[side]


class ConfigureWeno(object):
    """ Object containing the parameters needed by the WENO reconstruction for a given side.

    :arg int k: WENO coefficient k, equal to scheme order = 2k - 1.
    :arg int side: Side of the WENO reconstruction, either -1 (left) or +1 (right)."""

    def __init__(self, k, side):
        self.side = side
        if self.side == -1:
            self.name, self.short_name, self.shift = 'left', 'L', 1
        elif self.side == 1:
            self.name, self.short_name, self.shift = 'right', 'R', 0
        self.k, self.side = k, side
        self.func_points = self.generate_left_right_points()
        # k passed explicitly as 2 sets of ENO coefficients are needed for the smoothness indicators
        self.c_rj = self.generate_eno_coefficients(k)
        self.c2_rj = self.generate_eno_coefficients(2*k-1)
        self.opt_weights = self.generate_optimal_weights()
        return

    @property
    def generate_symbolic_function_points(self):
        """ Creates a set of symbolic function points f[0], f[1], ... to be used in the
        WENO reconstruction, and a dictionary linking each function point to its index.

        :returns symbolic_functions: The symbolic function points.
        :returns symbolic_function_dictionary: A dictionary of index:symbolic function pairs."""
        f = IndexedBase('f')
        symbolic_functions = []
        symbolic_function_dictionary = {}
        for p in set(self.func_points):
            symbolic_functions.append(f[p])
            symbolic_function_dictionary[p] = symbolic_functions[-1]
        return symbolic_functions, symbolic_function_dictionary

    def generate_eno_coefficients(self, k):
        """ Generates the c_rj ENO reconstruction coefficients given in Table 2.1
        of 'Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes
        for Hyperbolic Conservation Laws' by Shu(1997). Computation is of equation (2.21).

        :arg int k: WENO coefficient k, equal to scheme order = 2k - 1.
        :returns: dict: c_rj: Dictionary in the form (r,j) : ENO coefficient c_rj."""
        if self.side == 1:
            d = 1
        elif self.side == -1:
            d = 0
        r_values, j_values = range(k), range(k)
        c_rj = {}
        for r in r_values:
            for j in j_values:
                c_rj_sum = 0
                for m in range(j+1, k+1):
                    top_sum = 0
                    bottom_product = 1
                    for l in [x for x in range(k+1) if x != m]:
                        top_product_total = 1
                        for q in [x for x in range(k+1) if (x != m and x != l)]:
                            top_product = r - q + d
                            top_product_total *= top_product
                        bottom_product *= m - l
                        top_sum += top_product_total
                    c_rj_sum += Rational(top_sum, bottom_product)
                c_rj[(r, j)] = c_rj_sum
        return c_rj

    def generate_optimal_weights(self):
        """ Creates the optimal weights d_r needed by the WENO scheme, as specified in Shu(1997).

        :returns: coeff_dict: A dictionary of the optimal weights. """
        # Store reconstruction coefficients to form the sum
        k, c_rj, c2_rj = self.k, self.c_rj, self.c2_rj
        opt_weights = [Symbol('d_%d' % r) for r in range(k)]
        r_values = [[k-1-j for j in range(m+1)][::-1] for m in range(k)]
        equations = []
        for i, indices in enumerate(r_values):
            variables = [opt_weights[r] for r in indices]
            coeffs = [c_rj[(r, j)] for j, r in enumerate(indices)]
            rhs = c2_rj[(k-1, i)]
            v = [variables[i]*coeffs[i] for i in range(len(indices))]
            equations += [sum(v) - rhs]
        # Solve the linear system to get the coefficients
        solution = solve(equations, opt_weights)
        coeff_dict = {}
        for i in range(k):
            coeff_dict[(0, i)] = solution[opt_weights[i]]
        return coeff_dict

    def generate_left_right_points(self):
        """ Evaluate the required function evaluation points for left and
        right WENO reconstructions.

        :returns: fn_points: A list of integers to index the function locations."""
        k = self.k
        if self.side == -1:
            c = 1
        elif self.side == 1:
            c = 0
        fn_points = []
        for r in range(k):
            fn_points += [-r+j+c for j in range(k)]
        return fn_points


class JS_smoothness(object):
    """ An object to hold the Jiang-Shu smoothness coefficients. Listed as Eq(3.1) in
    'WENO coefficient k, equal to scheme order = 2k - 1.' Shu(1996)."""

    def __init__(self, k):
        self.k = k
        self.smooth_coeffs = self.generate_smoothness_coefficients()
        return

    def generate_smoothness_coefficients(self):
        """Extracts the JS smoothness coefficients.

        :returns: Dictionary of coefficients used in the generation of the smoothness indicators
        calculated in the 'generate_function_smoothness_indicators' routine."""
        k = self.k
        smooth_coeffs = {}
        z, dz = Symbol('z'), Symbol('dz')
        for r in range(k):
            # Generate z, y values to create interpolating polynomial
            dz_values = [dz*i for i in range(-r, k+1-r)]
            nodes = [Symbol('fn[i%+d]' % i) for i in range(-r, k-r)]
            funcs = [0]
            for i in range(k):
                funcs.append(funcs[-1] + dz*nodes[i])
            # Perform Lagrange interpolation
            lagrange_interp = interpolating_poly(k+1, z, dz_values, funcs).diff(z)
            # Loop over the derivatives of the interpolating polynomial
            total = 0
            for l in range(1, k):
                q = (lagrange_interp.diff(z, l))**2
                # Perform the integration and multiply by h^(2*l-1) over cell
                q = integrate(q.as_poly(z), z) * dz**(2*l-1)
                total += (q.subs(z, dz) - q.subs(z, 0))
            done = []
            # Save the coefficients of the smoothness indicator
            for m in range(0, 2*k-1):
                for n in range(0, 2*k-1):
                    func_product = Symbol('fn[i%+d]' % (-r+m))*Symbol('fn[i%+d]' % (-r+n))
                    if func_product not in done:
                        c = total.coeff(func_product)
                    if c != 0:
                        smooth_coeffs[(r, m, n)] = c
                    done.append(func_product)
        return smooth_coeffs

    def generate_function_smoothness_indicators(self, fn):
        """ Generates the standard Jiang-Shu smoothness indicators as polynomials in the
        symbolic placeholder functions 'f'.

        :arg object fn: The reconstruction variable object to contain the smoothness indicators."""
        if isinstance(fn, LeftWenoReconstructionVariable):
            shift = 1
        elif isinstance(fn, RightWenoReconstructionVariable):
            shift = 0
        k = self.k
        # Compute the smoothness indicator and alpha
        for r in range(k):
            fn.smoothness_symbols += [Symbol('beta_%d' % r)]
            local_smoothness = 0
            for m in range(0, 2*k-1):
                for n in range(0, 2*k-1):
                    beta = self.smooth_coeffs.get((r, m, n))
                    if beta is not None:
                        shift_indices = [-r+m+shift, -r+n+shift]
                        func_product = fn.function_stencil_dictionary[shift_indices[0]]*fn.function_stencil_dictionary[shift_indices[1]]
                        local_smoothness += beta*func_product
            fn.smoothness_indicators += [horner(local_smoothness)]
        return


class WenoReconstructionVariable(object):
    """ Reconstruction variable object to hold the quantities required for WENO.

    :arg str name: Name of the reconstruction, either left or right."""

    def __init__(self, name):
        self.name = name
        self.smoothness_indicators = []
        self.smoothness_symbols = []
        self.alpha_evaluated = []
        self.alpha_symbols = []
        self.inv_alpha_sum_symbols = []
        self.inv_alpha_sum_evaluated = []
        self.omega_evaluated = []
        self.omega_symbols = []
        self.function_stencil_dictionary = {}
        self.reconstructed_expression = None
        self.reconstructed_symbol = GridVariable('%s' % (name))
        return

    def update_quantities(self, original):
        """ Updates the quantities required by WENO in the reconstruction variable.

        :arg object original: Reconstruction object variable, either left or right reconstruction."""
        self.smoothness_symbols += [GridVariable('%s' % (s)) for s in original.smoothness_symbols]
        self.alpha_symbols += [GridVariable('%s' % (s)) for s in original.alpha_symbols]
        self.inv_alpha_sum_symbols += [GridVariable('%s' % (s)) for s in original.inv_alpha_sum_symbols]
        self.omega_symbols += [GridVariable('%s' % (s)) for s in original.omega_symbols]

        subs_dict = dict(zip(original.smoothness_symbols+original.alpha_symbols+original.inv_alpha_sum_symbols + original.omega_symbols, self.smoothness_symbols+self.alpha_symbols+self.inv_alpha_sum_symbols+self.omega_symbols))
        for key, value in original.function_stencil_dictionary.iteritems():
            subs_dict[value] = self.function_stencil_dictionary[key]

        self.smoothness_indicators = [s.subs(subs_dict) for s in original.smoothness_indicators]
        self.alpha_evaluated = [s.subs(subs_dict) for s in original.alpha_evaluated]
        self.inv_alpha_sum_evaluated = [s.subs(subs_dict) for s in original.inv_alpha_sum_evaluated]
        self.omega_evaluated = [s.subs(subs_dict) for s in original.omega_evaluated]
        self.reconstructed_expression = original.reconstructed_expression.subs(subs_dict)
        return

    def evaluate_quantities(self):
        all_symbols = self.smoothness_symbols + self.alpha_symbols + self.inv_alpha_sum_symbols + self.omega_symbols
        all_evaluations = self.smoothness_indicators + self.alpha_evaluated + self.inv_alpha_sum_evaluated + self.omega_evaluated
        final_equations = []
        for no, value in enumerate(all_symbols):
            final_equations += [OpenSBLIEq(value, all_evaluations[no])]
        self.final_equations = final_equations
        rv = self.reconstructed_symbol
        if self.settings.has_key("combine_reconstructions") and self.settings["combine_reconstructions"]:
            self.final_equations += [OpenSBLIEq(rv, rv + self.reconstructed_expression)]
        else:
            self.final_equations += [OpenSBLIEq(rv, self.reconstructed_expression)]
        return


class LeftWenoReconstructionVariable(WenoReconstructionVariable):
    """ Reconstruction object for the left reconstruction.

    :arg str name: 'left' """

    def __init__(self, name):
        WenoReconstructionVariable.__init__(self, name)
        return


class RightWenoReconstructionVariable(WenoReconstructionVariable):
    """ Reconstruction object for the right reconstruction.

    :arg str name: 'right' """

    def __init__(self, name):
        WenoReconstructionVariable.__init__(self, name)
        return


class WenoZ(object):
    def __init__(self, k):
        self.k = k
        self.eps = 1e-16
        return

    def global_smoothness_indicator(self, RV):
        """ Creates the WENO-Z global smoothness indicator.
        Formulae taken from 'Accuracy of the weighted essentially
        non-oscillatory conservative finite difference schemes.
        W-S. Don, R. Borges (2013). http://dx.doi.org/10.1016/j.jcp.2013.05.018.'

        :arg object RV: The reconstruction variable object."""
        b = RV.smoothness_symbols
        k = self.k  # WENO coefficient, order of scheme = 2k-1
        if k == 2:
            tau = Abs(b[0]-b[1])
        elif k == 3:
            tau = Abs(b[0]-b[2])
        elif k == 4:
            tau = Abs(b[0]+3*b[1]-3*b[2]-b[3])
        elif k == 5:
            tau = Abs(b[0]+2*b[1]-6*b[2]+2*b[3]+b[4])
        elif k == 6:
            tau = Abs(b[0]+b[1]-8*b[2]+8*b[3]-b[4]-b[5])
        elif k == 7:
            tau = Abs(b[0]+36*b[1]+135*b[2]-135*b[4]-36*b[5]-b[6])
        else:
            raise ValueError("WENO-Z global smoothness indicators only implemented up to WENO order 13.")
        return tau

    def generate_alphas(self, RV, WenoConfig):
        """ Create the alpha terms for the non-linear WENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        tau = self.global_smoothness_indicator(RV)
        p = 2
        for r in range(self.k):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append(WenoConfig.opt_weights.get((0, r))*(Float(1) + (tau/(self.eps + RV.smoothness_symbols[r]))**p))
        return

    def generate_omegas(self, RV, WenoConfig):
        """ Create the omega terms for the non-linear WENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        inv_alpha_sum_symbols = [Symbol('inv_alpha_sum')]
        inv_alpha_sum_evaluated = [S.One/sum(RV.alpha_symbols)]
        for r in range(self.k):
            RV.omega_symbols += [Symbol('omega_%d' % r)]
            RV.omega_evaluated += [RV.alpha_symbols[r]*inv_alpha_sum_symbols[0]]
        RV.inv_alpha_sum_symbols = inv_alpha_sum_symbols
        RV.inv_alpha_sum_evaluated = inv_alpha_sum_evaluated
        return


class WenoJS(object):
    def __init__(self, k):
        self.k = k
        self.eps = 1e-6
        return

    def generate_alphas(self, RV, WenoConfig):
        """ Create the alpha terms for the non-linear WENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        for r in range(self.k):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append(WenoConfig.opt_weights.get((0, r))/(self.eps+RV.smoothness_symbols[r])**2)
        return

    def generate_omegas(self, RV, WenoConfig):
        """ Create the omega terms for the non-linear WENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        inv_alpha_sum_symbols = [Symbol('inv_alpha_sum')]
        inv_alpha_sum_evaluated = [S.One/sum(RV.alpha_symbols)]
        for r in range(self.k):
            RV.omega_symbols += [Symbol('omega_%d' % r)]
            RV.omega_evaluated += [RV.alpha_symbols[r]*inv_alpha_sum_symbols[0]]
        RV.inv_alpha_sum_symbols = inv_alpha_sum_symbols
        RV.inv_alpha_sum_evaluated = inv_alpha_sum_evaluated
        return


class Weno(Scheme, ShockCapturing):
    """ Main WENO class. Performs the Jiang-Shu WENO reconstruction procedure. Refer to the reference:
    'Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory schemes for Hyperbolic Conservation
    laws. Shu (1997). The alpha & omega quantities follow the description from this paper.

    :arg int order: Numerical order of the WENO scheme (3,5,...). """

    def __init__(self, order, formulation):
        Scheme.__init__(self, "WenoDerivative", order)
        self.schemetype = "Spatial"
        if not isinstance(order, int):
            raise TypeError("Weno order should be an integer, if using WenoZ: pass formulation ='Z'")
        self.k = int(0.5*(order+1))
        self.order = order
        if formulation.upper() == 'Z':
            WT = WenoZ(self.k)
            print("A WENO-Z scheme of order %s is being used for shock capturing." % str(self.order))
        elif formulation.upper() == 'JS':  # Default to WENO-JS if no WENO scheme type provided
            WT = WenoJS(self.k)
            print("A WENO-JS scheme of order %s is being used for shock capturing." % str(self.order))
        else:
            raise NotImplementedError("Only WENO-Z and WENO-JS schemes are implemented.")
        JS = JS_smoothness(self.k)
        self.halotype = WenoHalos(self.order)
        self.required_constituent_relations_symbols = {}
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        self.reconstruction_classes = [RightWenoReconstructionVariable('right'), LeftWenoReconstructionVariable('left')]
        # Populate the quantities required by WENO for the left and right reconstruction variable.
        for no, side in enumerate([1, -1]):
            WenoConfig = ConfigureWeno(self.k, side)
            RV = self.reconstruction_classes[no]
            RV.func_points = sorted(set(WenoConfig.func_points))
            RV.stencil_points, RV.function_stencil_dictionary = WenoConfig.generate_symbolic_function_points
            JS.generate_function_smoothness_indicators(RV)
            WT.generate_alphas(RV, WenoConfig)
            WT.generate_omegas(RV, WenoConfig)
            # Final reconstruction is the same for WENO-JS and WENO-Z
            self.generate_reconstruction(RV, WenoConfig)
            self.reconstruction_classes[no] = RV
        return

    def reconstruction_halotype(self, order, reconstruction=True):
        return WenoHalos(order, reconstruction)

    def group_by_direction(self, eqs):
        """ Groups the input equations by the direction (x0, x1, ...) they depend upon.

        :arg list eqs: List of equations to group by direction.
        :returns: dict: grouped: Dictionary of {direction: equations} key, value pairs for equations grouped by direction."""
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(WenoDerivative))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped

    def generate_reconstruction(self, RV, WenoConfig):
        """ Create the final WENO stencil by summing the stencil points, ENO coefficients and WENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object WenoConfig: Configuration settings for a reconstruction of either left or right."""
        stencil = 0
        k, c_rj = self.k, WenoConfig.c_rj
        for r in range(k):
            eno_expansion = [c_rj.get((r, j))*RV.function_stencil_dictionary[WenoConfig.func_points[r*k+j]] for j in range(k)]
            stencil += RV.omega_symbols[r]*sum(eno_expansion)
        RV.reconstructed_expression = stencil
        return


class LLFWeno(LLFCharacteristic, Weno):
    """ Performs the Local Lax Friedrichs flux splitting with a WENO scheme.

    :arg int order: Order of the WENO/TENO scheme.
    :arg object physics: Physics object, defaults to NSPhysics.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging. """

    def __init__(self, order, physics=None, averaging=None, formulation="JS"):
        print("Local Lax-Friedrich flux splitting.")
        LLFCharacteristic.__init__(self, physics, averaging)
        Weno.__init__(self, order, formulation)
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
                # Get the pre, interpolations and post equations for characteristic reconstruction
                pre_process, interpolated, post_process = self.get_characteristic_equations(direction, derivatives, solution_vector, block)
                # Add the equations to the kernel and add the kernel to SimulationEquations
                kernel.add_equation(pre_process + interpolated + post_process)
                type_of_eq.Kernels += [kernel]
            # Generate kernels for the constituent relations
            if grouped:
                constituent_relations = self.generate_constituent_relations_kernels(block)
                type_of_eq.Kernels += [self.evaluate_residuals(block, eqs, all_derivatives_evaluated_locally)]
                constituent_relations = self.check_constituent_relations(block, eqs, constituent_relations)
            return constituent_relations


class RFWeno(RFCharacteristic, Weno):
    """ Performs the Local Lax Friedrichs flux splitting with a WENO scheme.

    :arg int order: Order of the WENO/TENO scheme.
    :arg object physics: Physics object, defaults to NSPhysics.
    :arg object averaging: The averaging procedure to be applied for characteristics, defaults to Simple averaging. """

    def __init__(self, order, physics=None, averaging=None, formulation="JS"):
        print("Roe-flux differencing.")
        RFCharacteristic.__init__(self, physics, averaging)
        Weno.__init__(self, order, formulation)
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
                # Get the pre, interpolations and post equations for characteristic reconstruction
                pre_process, interpolated, post_process = self.get_characteristic_equations(direction, derivatives, solution_vector, block)
                # Add the equations to the kernel and add the kernel to SimulationEquations
                kernel.add_equation(pre_process + interpolated + post_process)
                type_of_eq.Kernels += [kernel]
            # Generate kernels for the constituent relations
            if grouped:
                constituent_relations = self.generate_constituent_relations_kernels(block)
                type_of_eq.Kernels += [self.evaluate_residuals(block, eqs, all_derivatives_evaluated_locally)]
                constituent_relations = self.check_constituent_relations(block, eqs, constituent_relations)
            return constituent_relations

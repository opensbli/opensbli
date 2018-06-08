from sympy import IndexedBase, Symbol, Rational, Abs, flatten, Max, horner, S, ceiling, pprint
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
            if self.order == 8:
                coeffs += [[Rational(-3, 12), Rational(13, 12), Rational(-23, 12), Rational(25, 12)],
                           [Rational(12, 60), Rational(77, 60), Rational(-43, 60), Rational(17, 60), Rational(-3, 60)]]
        elif side == -1:
            coeffs = [[Rational(2, 6), Rational(5, 6), Rational(-1, 6)], [Rational(-1, 6), Rational(5, 6), Rational(2, 6)],
                      [Rational(11, 6), Rational(-7, 6), Rational(1, 3)]]
            if self.order in [6, 8]:
                coeffs += [[Rational(1, 12), Rational(-5, 12), Rational(13, 12), Rational(3, 12)]]
            if self.order == 8:
                coeffs += [[Rational(25, 12), Rational(-23, 12), Rational(13, 12), Rational(-1, 4)],
                           [Rational(-1, 20), Rational(17, 60), Rational(-43, 60), Rational(77, 60), Rational(12, 60)]]
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
            elif order == 8:
                opt_coeffs = [0.4336570089737348, 0.2193140179474722, 0.07144766367542149, 0.1302093452983125,
                              0.03089532735084351, 0.1144766367542177]
        elif not self.optimized:
            if order == 5:
                opt_coeffs = [Rational(6, 10), Rational(3, 10), Rational(1, 10)]
            elif order == 6:
                opt_coeffs = [Rational(9, 20), Rational(6, 20), Rational(1, 20), Rational(4, 20)]
            elif order == 8:
                opt_coeffs = [Rational(30, 70), Rational(18, 70), Rational(4, 70), Rational(12, 70), Rational(1, 70), Rational(5, 70)]
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
            # Add Stencils [4, 5] for TENO8
            if self.order == 8:
                smoothness_indicators += [Rational(1, 36)*(-2*fns[-3]+9*fns[-2]-18*fns[-1]+11*fns[0])**2 +
                                          Rational(13, 12)*(-1*fns[-3]+4*fns[-2]-5*fns[-1]+2*fns[0])**2 +
                                          Rational(781, 720)*(-1*fns[-3]+3*fns[-2]-3*fns[-1]+fns[0])**2]
                # Use horner on this one, reduces number of operations
                smoothness_indicators += [horner(Rational(1, 144)*(-25*fns[0]+48*fns[1]-36*fns[2]+16*fns[3]-3*fns[4])**2 +
                                                 Rational(13, 1728)*(35*fns[0]-104*fns[1]+114*fns[2]-56*fns[3]+11*fns[4])**2 +
                                                 Rational(781, 2880)*(-5*fns[0]+18*fns[1]-24*fns[2]+14*fns[3]-3*fns[4])**2 -
                                                 Rational(1, 4320)*(35*fns[0]-104*fns[1]+114*fns[2]-56*fns[3]+11*fns[4]) *
                                                 (fns[0]-4*fns[1]+6*fns[2]-4*fns[3]+fns[4]) + Rational(32803, 30240) *
                                                 (fns[0]-4*fns[1]+6*fns[2]-4*fns[3]+fns[4])**2)]
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
            if self.order == 8:
                # [0, 1, 2, 3]
                # smoothness_indicators += [Rational(1,36)*(-2*fns[0]-3*fns[1]+6*fns[2]-fns[3])**2 + Rational(13,12)*(fns[0]-2*fns[1]+fns[2])**2 +\
                #                          + Rational(1043,960)*(-fns[0]+3*fns[1]-3*fns[2]+fns[3])**2 + Rational(1,432)*(-2*fns[0]-3*fns[1]+6*fns[2]-fns[3])*\
                #                          (-fns[0]+3*fns[1]-3*fns[2]+fns[3])]
                # smoothness_indicators += [horner(Rational(547, 240)*fns[0]**2 - Rational(1261, 120)*fns[0]*fns[1]+Rational(961, 120)*fns[0]*fns[2] - Rational(247, 120)*fns[0]*fns[3]
                #                                  + Rational(3443, 240)*fns[1]**2 - Rational(2983, 120)*fns[1]*fns[2] + Rational(267, 40)*fns[1]*fns[3]+Rational(2843, 240)*fns[2]**2
                #                                  - Rational(821, 120)*fns[2]*fns[3] + Rational(89, 80)*fns[3]**2)]
                smoothness_indicators += [(7043*fns[3]/240 - 647*fns[4]/40)*fns[3] + (11003*fns[2]/240 - 8623*fns[3]/120 + 2321*fns[4]/120)*fns[2] +
                                          (2107*fns[1]/240 - 1567*fns[2]/40 + 3521*fns[3]/120 - 309*fns[4]/40)*fns[1] + 547*fns[4]**2/240]
                # Factored , [-3, -2, -1, 0, 1]
                smoothness_indicators += [horner(Rational(11329, 2520)*fns[-3]**2-Rational(208501, 5040)*fns[-3]*fns[-2]+Rational(121621, 1680)*fns[-3]*fns[-1]
                                                 - Rational(288007, 5040)*fns[-3]*fns[0]+Rational(86329, 5040)*fns[-3]*fns[1]+Rational(482963, 5040)*fns[-2]**2
                                                 - Rational(142033, 420)*fns[-2]*fns[-1]+Rational(679229, 2520)*fns[-2]*fns[0]-Rational(411487, 5040)*fns[-2]*fns[1]
                                                 + Rational(507131, 1680)*fns[-1]**2 - Rational(68391, 140)*fns[-1]*fns[0]+Rational(252941, 1680)*fns[-1]*fns[1]
                                                 + Rational(1020563, 5040)*fns[0]**2 - Rational(649501, 5040)*fns[0]*fns[1]+Rational(53959, 2520)*fns[1]**2)]
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


class Teno8(object):
    """ Base class for the 8th order TENO scheme, still work in progress."""

    def __init__(self):
        return

    def global_smoothness_indicator(self, RV):
        # Global smoothness indicator used in tau_8 for TENO8
        points = [-3, -2, -1, 0, 1, 2, 3, 4]
        f = IndexedBase('f')
        symbolic_functions = []
        fns = {}
        for p in points:
            symbolic_functions.append(f[p])
            fns[p] = symbolic_functions[-1]
        # tau_8 = fns[4]*(75349098471*fns[4]-1078504915264*fns[3]+3263178215782*fns[2]
        #                 - 5401061230160*fns[1]+5274436892970*fns[0]-3038037798592*fns[-1]+956371298594*fns[-2]-127080660272*fns[-3])+fns[3]*(3944861897609*fns[3]-24347015748304*fns[2]
        #                                                                                                                                      + 41008808432890*fns[1]-40666174667520*fns[0]+23740865961334*fns[-1]-7563868580208*fns[-2]
        #                                                                                                                                      + 1016165721854*fns[-3])+fns[2]*(38329064547231*fns[2]-131672853704480*fns[1]+132979856899250*fns[0]-78915800051952*fns[-1]+25505661974314*fns[-2] - 3471156679072*fns[-3])\
        #     + fns[1]*(115451981835025*fns[1]-238079153652400*fns[0]+144094750348910*fns[-1]
        #               - 47407534412640*fns[-2]+6553080547830*fns[-3])+fns[0]*(125494539510175*fns[0]
        #                                                                       - 155373333547520*fns[-1]+52241614797670*fns[-2]-7366325742800*fns[-3]) \
        #     + fns[-1]*(49287325751121*fns[-1]-33999931981264*fns[-2]+4916835566842*fns[-3]) \
        #     + fns[-2]*(6033767706599*fns[-2]-1799848509664*fns[-3])+139164877641*fns[-3]*fns[-3]
        # tau_8 = Rational(1, 62270208000)*(tau_8)
        tau_8 = [(3944861897609*fns[3]/62270208000 - 2407377043*fns[4]/138996000)*fns[3] + (12780967457077*fns[2]/20756736000 - 1521688484269*fns[3]/3891888000 +
                                                                                            1631589107891*fns[4]/31135104000)*fns[2] + (420341161931*fns[1]/226437120 - 74851467823*fns[2]/35380800 + 4100880843289*fns[3]/6227020800 - 67513265377*fns[4]/778377600)*fns[1] +
                 (457249528517*fns[0]/226437120 - 595915721251*fns[1]/155675520 + 532071643661*fns[2]/249080832 - 10590149653*fns[3]/16216200 + 25116366157*fns[4]/296524800)*fns[0] +
                 (46388292547*fns[-3]/20756736000 - 18415814357*fns[0]/155675520 + 72812006087*fns[1]/691891200 - 108473646221*fns[2]/1945944000 + 508082860927*fns[3]/31135104000 -
                  7942541267*fns[4]/3891888000)*fns[-3] + (6047605530599*fns[-2]/62270208000 - 56245265927*fns[-3]/1945944000 + 5227966881367*fns[0]/6227020800 - 98765696693*fns[1]/129729600 +
                                                           12752830987157*fns[2]/31135104000 - 157580595421*fns[3]/1297296000 + 478185649297*fns[4]/31135104000)*fns[-2] + (16476387815707*fns[-1]/20756736000 - 2129103852829*fns[-2]/3891888000 +
                                                                                                                                                                            2458417783421*fns[-3]/31135104000 - 5527715497*fns[0]/2211300 + 14416393946891*fns[1]/6227020800 - 1644079167749*fns[2]/1297296000 + 11870432980667*fns[3]/31135104000 -
                                                                                                                                                                            47469340603*fns[4]/972972000)*fns[-1] + 25116366157*fns[4]**2/20756736000]
        RV.tau_8_symbol = [Symbol('tau')]
        RV.tau_8_evaluated = tau_8
        return

    def generate_alphas(self, RV, TC):
        """ Create the alpha terms for the non-linear TENO weights.

        :arg object RV: The reconstruction variable object.
        :arg object TC: Configuration settings for a reconstruction of either left or right."""
        # Scale separation parameters
        C, q = S.One, 6
        RV.tau_8_symbol = [Symbol('tau')]
        tau_8 = RV.tau_8_symbol[0] - Rational(1, 6)*(RV.smoothness_symbols[1] + RV.smoothness_symbols[2] + 4*RV.smoothness_symbols[0])
        for r in range(TC.n_stencils):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append((C + (Abs(tau_8)/(self.eps + RV.smoothness_symbols[r])))**q)
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
        elif order == 8:
            WT = Teno8()
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

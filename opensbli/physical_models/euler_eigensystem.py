from opensbli.core.opensbliobjects import CoordinateObject
from sympy import diag, sqrt, Eq, eye, Rational, pprint
from opensbli.core.opensbliobjects import ConstantObject, EinsteinTerm
from sympy.parsing.sympy_parser import parse_expr
# from opensbli.core.metrics import MetricsEquation


class EulerEquations(object):
    def __init__(self, ndim, **kwargs):
        self.ndim = ndim
        return

    def apply_direction(self, direction):
        # Dictionaries to store ev, REV and LEV for each direction
        ev_dict, LEV_dict, REV_dict = {}, {}, {}
        # Metric symbols from block
        met_symbols = self.met_symbols
        # Metric terms for this direction to substitute into the matrix
        terms = [EinsteinTerm('k%d' % i) for i in range(self.ndim)]
        metric_values = [met_symbols[direction,i] for i in range(self.ndim)]
        subs_dict = dict([(x,y) for (x,y) in zip(terms, metric_values)])
        # Scaling factor based on metrics
        factor = sum([met_symbols[direction, i]**2 for i in range(self.ndim)])**(Rational(1,2))
        subs_dict[EinsteinTerm('k')] = factor
        def g(x):
            return x.subs(subs_dict, evaluate=False)
        ev_dict[direction] = diag(*list(self.ev.applyfunc(g)))
        LEV_dict[direction] = self.LEV.applyfunc(g)
        REV_dict[direction] = self.REV.applyfunc(g)
        return ev_dict, LEV_dict, REV_dict

    def generate_eig_system(self, block):
        ndim = self.ndim
        # Check if block has metrics equations:
        self.met_symbols = block.fd_metrics
        local_dict = {'Symbol': EinsteinTerm, 'gama': ConstantObject('gama')}

        if ndim == 1:
            matrix_symbols = ['H']
            matrix_formulae = ['a**2/(gama-1) + u0**2/2']
            matrix_symbols = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae]
            ev = 'diag([u0-a, u0, u0+a])'
            REV = 'Matrix([[1,1,1], [u0-a,u0,u0+a], [H-u0*a,u0**2 /2,H+u0*a]])'
            LEV = 'Matrix([[ u0*(2*H + a*u0 - u0**2)/(2*a*(2*H - u0**2)), (-H - a*u0 + u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)],[2*(H - u0**2)/(2*H - u0**2),2*u0/(2*H - u0**2), 2/(-2*H + u0**2)],[u0*(-2*H + a*u0 + u0**2)/(2*a*(2*H - u0**2)),  (H - a*u0 - u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)]])'
            ev = parse_expr(ev, local_dict=local_dict, evaluate=False)
            REV = parse_expr(REV, local_dict=local_dict, evaluate=False)
            LEV = parse_expr(LEV, local_dict=local_dict, evaluate=False)

            subs_dict = dict(zip(matrix_symbols, matrix_formulae))

            def f(x):
                return x.subs(subs_dict)
            ev = ev.applyfunc(f)
            REV = REV.applyfunc(f)
            LEV = LEV.applyfunc(f)
            ev_dict = {}
            REV_dict = {}
            LEV_dict = {}
            definitions = {}

            subs_list = [{EinsteinTerm('un_0'): 1, EinsteinTerm('un_1'): 0},
                         {EinsteinTerm('un_0'): 0, EinsteinTerm('un_1'): 1}]
            equations = [Eq(a, b) for a, b in subs_dict.items()]
            eq_directions = {}
            for no, coordinate in enumerate(coordinates):
                direction = coordinate.direction

                def g(x):
                    return x.subs(subs_list[no]).simplify()
                definitions[direction] = ev.applyfunc(g).atoms(EinsteinTerm).union(REV.applyfunc(g).atoms(EinsteinTerm)).union(LEV.applyfunc(g).atoms(EinsteinTerm))
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)
                eq_directions[direction] = [e.subs(subs_list[no]) for e in equations]

        if ndim == 2:
            matrix_symbols = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            matrix_formulae = ['rho/(a*sqrt(2.0))', '1/(rho*a*sqrt(2.0))', 'k0*u0 + k1*u1', '(gama-1)*((u0**2 + u1**2)/2)', 'sqrt(k0**2 + k1**2)']
            matrix_symbols_2 = ['k0_t', 'k1_t', 'theta_t', 'U']
            matrix_formulae_2 = ['k0/k', 'k1/k', 'k0_t*u0 + k1_t*u1', 'k0*u0+k1*u1']

            matrix_symbols = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae]
            matrix_symbols_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols_2]
            matrix_formulae_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae_2]
            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            subs_dict2 = dict(zip(matrix_symbols_2, matrix_formulae_2))

            ev = 'diag([U, U, U+a*k, U-a*k])'
            REV = 'Matrix([[1,0,alpha,alpha], [u0,k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [u1,-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [phi_sq/(gama-1),rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + a*theta_t),alpha*((phi_sq+a**2)/(gama-1) - a*theta_t)]])'
            LEV = 'Matrix([[1-phi_sq/a**2,(gama-1)*u0/a**2,(gama-1)*u1/a**2,-(gama-1)/a**2], [-(k1_t*u0-k0_t*u1)/rho,k1_t/rho,-k0_t/rho,0], [bta*(phi_sq-a*theta_t),bta*(k0_t*a-(gama-1)*u0),bta*(k1_t*a-(gama-1)*u1),bta*(gama-1)], [bta*(phi_sq+a*theta_t),-bta*(k0_t*a+(gama-1)*u0),-bta*(k1_t*a+(gama-1)*u1),bta*(gama-1)]])'
            ev = parse_expr(ev, local_dict=local_dict, evaluate=False)
            REV = parse_expr(REV, local_dict=local_dict, evaluate=False)
            LEV = parse_expr(LEV, local_dict=local_dict, evaluate=False)
            # Apply the sub
            def f(x):
                return x.subs(subs_dict, evaluate=False)
            def g(x):
                return x.subs(subs_dict2, evaluate=False)
            self.ev = ev.applyfunc(f).applyfunc(g).applyfunc(g)
            self.REV = REV.applyfunc(f).applyfunc(g).applyfunc(g)
            self.LEV = LEV.applyfunc(f).applyfunc(g).applyfunc(g)
            return

        elif ndim == 3:
            matrix_symbols = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            matrix_formulae = ['rho/(a*sqrt(2.0))', '1/(rho*a*sqrt(2.0))', 'k0*u0 + k1*u1 + k2*u2', '(gama-1)*((u0**2 + u1**2 + u2**2)/2)', 'sqrt(k0**2 + k1**2 + k2**2)']
            matrix_symbols_2 = ['k0_t', 'k1_t', 'k2_t', 'theta_t', 'U']
            matrix_formulae_2 = ['k0/k', 'k1/k', 'k2/k', 'k0_t*u0 + k1_t*u1 + k2_t*u2', 'k0*u0+k1*u1+k2*u2']

            matrix_symbols = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae]
            matrix_symbols_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols_2]
            matrix_formulae_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae_2]
            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            subs_dict2 = dict(zip(matrix_symbols_2, matrix_formulae_2))

            ev = 'diag([U, U, U, U+a*k, U-a*k])'
            REV = 'Matrix([[k0_t,k1_t,k2_t,alpha,alpha], [k0_t*u0,k1_t*u0-k2_t*rho,k2_t*u0+k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [k0_t*u1+k2_t*rho,k1_t*u1,k2_t*u1-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [k0_t*u2-k1_t*rho,k1_t*u2+k0_t*rho,k2_t*u2,alpha*(u2+k2_t*a),alpha*(u2-k2_t*a)], [k0_t*phi_sq/(gama-1) + rho*(k2_t*u1-k1_t*u2),k1_t*phi_sq/(gama-1) + rho*(k0_t*u2-k2_t*u0),k2_t*phi_sq/(gama-1) + rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + theta_t*a),alpha*((phi_sq+a**2)/(gama-1) - theta_t*a)]])'
            LEV = 'Matrix([[k0_t*(1-phi_sq/a**2)-(k2_t*u1-k1_t*u2)/rho,k0_t*(gama-1)*u0/a**2,k0_t*(gama-1)*u1/a**2 + k2_t/rho,k0_t*(gama-1)*u2/a**2 - k1_t/rho,-k0_t*(gama-1)/a**2], [k1_t*(1-phi_sq/a**2)-(k0_t*u2-k2_t*u0)/rho,k1_t*(gama-1)*u0/a**2 - k2_t/rho,k1_t*(gama-1)*u1/a**2,k1_t*(gama-1)*u2/a**2 + k0_t/rho,-k1_t*(gama-1)/a**2], [k2_t*(1-phi_sq/a**2) - (k1_t*u0-k0_t*u1)/rho,k2_t*(gama-1)*u0/a**2 + k1_t/rho,k2_t*(gama-1)*u1/a**2 - k0_t/rho,k2_t*(gama-1)*u2/a**2,-k2_t*(gama-1)/a**2], [bta*(phi_sq-theta_t*a),-bta*((gama-1)*u0-k0_t*a),-bta*((gama-1)*u1-k1_t*a),-bta*((gama-1)*u2 - k2_t*a),bta*(gama-1)], [bta*(phi_sq+theta_t*a),-bta*((gama-1)*u0+k0_t*a),-bta*((gama-1)*u1+k1_t*a),-bta*((gama-1)*u2+k2_t*a),bta*(gama-1)]])'
            ev = parse_expr(ev, local_dict=local_dict, evaluate=False)
            REV = parse_expr(REV, local_dict=local_dict, evaluate=False)
            LEV = parse_expr(LEV, local_dict=local_dict, evaluate=False)

            # Apply the sub
            def f(x):
                return x.subs(subs_dict, evaluate=False)

            def g(x):
                return x.subs(subs_dict2, evaluate=False)
            self.ev = ev.applyfunc(f).applyfunc(g).applyfunc(g)
            self.REV = REV.applyfunc(f).applyfunc(g).applyfunc(g)
            self.LEV = LEV.applyfunc(f).applyfunc(g).applyfunc(g)
        return

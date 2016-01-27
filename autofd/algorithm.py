''' This contains the routines used for algorithm generation
'''
from sympy import *
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_application)
transformations = standard_transformations + (implicit_application,)
import os

# AutoFD functions
from .utils import *
from .codegen_utils import *
from .fortran_subroutine import *

BUILD_DIR = os.getcwd()

import logging
LOG = logging.getLogger(__name__)


class AlgorithmInput(object):

    def __init__(self):
        self.time_sch = []  # Temporal discretisation scheme (e.g. RK3).
        self.time_ooa = []  # Temporal order of accuracy
        self.const_timestep = True  # Constant timestep.
        self.sp_sch = []  # Spatial discretisation scheme
        self.sp_ooa = []  # Spatial order of accuracy
        self.wrkarr = []  # Work array
        self.ceval = []  # Evaluation count
        self.bcs = []  # Boundary conditions
        self.lang = []  # Generated code's language
        self.multiblock = False  # Multi-block
        self.nblocks = 1  # Number of blocks
        self.expfilter = []  # Explicit filter
        return


def read_alg(read_file):
    comm_lineno = []
    # Get all the comments in the file
    for ind, line in enumerate(read_file):
        if line[0] == '#':
            comm_lineno.append(ind)
    alg_inp = AlgorithmInput()

    temp = read_file[comm_lineno[0]+1:comm_lineno[1]]
    if len(temp) > 1:
        raise ValueError('Temporal scheme is defined wrongly in algorithm text file')
    temp = temp[0].split(',')
    if temp[0] != 'RK':
        raise ValueError('Implement %s time stepping scheme' % temp[0])
    alg_inp.time_sch = temp[0]
    alg_inp.time_ooa = int(temp[1])
    temp = read_file[comm_lineno[1]+1:comm_lineno[2]]
    if len(temp) > 1:
        raise ValueError('Time step can be const or local only')
    if temp[0] == 'const':
        alg_inp.const_timestep = True
    elif temp[0] == 'local':
        alg_inp.const_timestep = False
    else:
        raise ValueError('Time step can be const or local only')
    temp = read_file[comm_lineno[2]+1:comm_lineno[3]]
    if len(temp) > 1:
        raise ValueError('Spatial scheme is defined wrongly in algorithm text file')
    temp = temp[0].split(',')
    if temp[0] != 'central_diff':
        raise ValueError('Implement %s spatial scheme' % temp[0])
    alg_inp.sp_sch = temp[0]
    alg_inp.sp_ooa = int(temp[1])
    temp = read_file[comm_lineno[3]+1:comm_lineno[4]]
    if len(temp) > 1:
        raise ValueError('ceval wrong')
    alg_inp.ceval = int(temp[0])
    temp = read_file[comm_lineno[4]+1:comm_lineno[5]]
    if len(temp) > 1:
        raise ValueError('ceval wrong')
    alg_inp.wrkarr = '%s%%d' % temp[0]
    temp = read_file[comm_lineno[5]+1:comm_lineno[6]]
    # Read the boundary conditions
    for te in temp:
        te = te.replace(' ', '')
        alg_inp.bcs = alg_inp.bcs + [tuple(te.strip().split(','))]
    # Language
    temp = read_file[comm_lineno[6]+1:comm_lineno[7]]
    alg_inp.lang = temp[0]
    # Multi-block stuff
    temp = read_file[comm_lineno[7]+1:comm_lineno[8]]
    if len(temp) > 1:
        raise ValueError('Blocks defined wrongly in the code')
    temp = temp[0]
    if temp == 'SB':
        alg_inp.multiblock = False
        alg_inp.nblocks = 1
    elif temp == 'MB':
        alg_inp.multiblock = True
        alg_inp.nblocks = Symbol('nblocks')
    else:
        raise ValueError('Blocks can be SB for single block or MB for Multi-block')

    # Spatial filtering of conservative varibles
    temp = read_file[comm_lineno[8]+1:comm_lineno[9]]
    temp = temp[0].split(',')
    for te in temp:
        te = te.strip()
        if te == 'T':
            alg_inp.expfilter.append(True)
        elif te == 'F':
            alg_inp.expfilter.append(False)
        else:
            raise ValueError('Filter can be T or F only')
    # Various checks performed on the algorithm

    return alg_inp


def perform_algorithm_checks(alg, inp):
    if len(alg.bcs) == inp.ndim:
        pass
    elif (len(alg.bcs) > inp.ndim):
        raise ValueError('There are more boundary conditions than the number of dimensions')
    elif (len(alg.bcs) < inp.ndim):
        raise ValueError('There are less boundary conditions than the number of dimensions')

    if len(alg.expfilter) == inp.ndim:
        pass
    elif (len(alg.expfilter) > inp.ndim):
        raise ValueError('There are more options for filter than number of dimension')
    elif (len(alg.expfilter) < inp.ndim):
        raise ValueError('There are less options for filter than number of dimension')

    return


def count_input_equations(variable, equations):
    """ Return the number of input equations containing a particular variable.
    
    :arg variable: the variable under consideration.
    :arg equations: the equations to search.
    :returns: the number of equations containing the variable.
    :rtype: int
    """

    count = 0
    for e in equations:
        count = count + e.count(variable)
    return count


def equations_to_dict(equations):
    """ Get the LHS and RHS of each equation, and return them in dictionary form.
    
    :arg equations: the equations to consider.
    :returns: a dictionary of (LHS, RHS) pairs.
    :rtype: dict
    """

    lhs = list(e.lhs for e in equations)
    rhs = list(e.rhs for e in equations)
    d = dict(zip(lhs, rhs))
    LOG.debug(d)
    return d


def sort_evals(inp, evald):
    inpcopy = [te for te in inp]
    out = []
    while inpcopy:
        # LOG.debug('in while loop')
        for te in inpcopy:
            ter = te.rhs.atoms(Indexed)
            if not ter.difference(set(evald)):
                out.append(te)
                inpcopy.remove(te)
                evald.append(te.lhs)
    return out, evald


def derform(inp, alg):
    ooa = alg.sp_ooa
    sch = alg.sp_sch
    points = []
    order = 2
    if sch == 'central_diff':
        points = list(i for i in range(-ooa/2, ooa/2+1))
        if len(points) < order+1:
            raise ValueError("Too few points for derivative of order %d" % order)
        wghts = finite_diff_weights(order, points, 0)
        inp.fd_wghts = wghts
        inp.fd_points = points
        for dim in range(inp.ndim):
            inp.halos.append(int(ooa/2))
            inp.gridhalo = inp.gridhalo + [inp.grid[dim+1]+2*inp.halos[dim]]
    else:
        raise ValueError('Implement %s scheme in derivative formulation' % sch)
    return


def apply_der(der, inp):
    order = len(der.args) - 1
    temp = []
    lwr = inp.block.lower
    upr = inp.block.upper
    for ind in inp.fd_points:
        outderivative = der.args[0]
        for d in der.args[0].atoms(Indexed):
            var = d.base
            index = d.indices
            dire = der.args[1]
            ind1 = str(index[0]).replace('%s' % str(dire), '%s' % str(dire + ind))
            ind1 = Idx(ind1, (lwr, upr))
            outderivative = outderivative.subs(d, var[ind1])
        temp = temp + [outderivative]
    derivative = 0
    N = len(inp.fd_points) - 1
    for nu in range(0, len(inp.fd_points)):
        derivative += inp.fd_wghts[order][N][nu]*temp[nu]
        # derivative = derivative.simplify()
        # derivative = ratsimp(derivative)
        # derivative = derivative.simplify()
        # derivative = ratsimp(derivative)
        # outderivative = outderivative.subs(d, derivative)
    if order == 1:
        derivative = derivative/(Symbol('d%s' % str(dire)))
        inp.const = inp.const + [Symbol('d%s' % str(dire))]
    elif order == 2:
        derivative = derivative/(Symbol('d%s' % str(dire))*Symbol('d%s' % str(dire)))
        inp.const = inp.const + [Symbol('d%s' % str(dire))]
    else:
        raise ValueError('Implement the order of derivative')
    # easydebug
    # derivative = derivative.simplify()
    # derivative = ratsimp(derivative)
    # end easydebug
    return derivative


class ExplicitFilter(object):

    """ Explicit filter implemented from the NASA PAPER. """

    def __init__(self):
        self.coeffs = {}
        self.coeffs[2] = [-1, 2, -1]
        self.coeffs[4] = [-1, 4, -6, 4, -1]

    def filter_equations(self, alg):
        """ Returns the filter equations. """
        filtereqs = []
        for dire, fil in enumerate(alg.expfilter):
            if fil:
                ooa = alg.sp_ooa
                alphaD = (-1)**(int(ooa/2) + 1) * (2) ** (-ooa)
                if self.coeffs.get(ooa):
                    muls = self.coeffs.get(ooa)
                    temp = list(con.base for con in inp.conser)
                    out = [con for con in inp.conser]
                    # out = 0
                    direction = Symbol('x%d' % dire)
                    for num, ind in enumerate(inp.fd_points):
                        repl = '%s' % str(direction + ind)
                        index = str(inp.varform).replace(str(direction), repl)
                        index = Idx(index, Symbol('nblock', integer=True))
                        for nv in range(len(temp)):
                            varib = temp[nv]
                            out[nv] += alphaD * muls[num] * varib[index]
                    for nv in range(len(temp)):
                        out[nv] = Eq(savedic.get(inp.conser[nv]), out[nv])
                    filtereqs.append(out)
            else:
                raise ValueError('Implement spatial filtering coefficients for order of accuracy %d', ooa)
        return filtereqs


class PreparedEquations(object):

    def __init__(self, eqs, forms, alginp):
        self.inpeq = (flatten(list(eq.expandedeq for eq in eqs)))
        self.forms = flatten(list(eq.expandedeq for eq in forms))
        var = flatten(list(eq.variables for eq in eqs))
        conser = flatten(list(eq.conser for eq in eqs))
        const = list(set(flatten(list(eq.const for eq in eqs))))
        totvar = list(set(conser + var))
        self.ndim = eqs[0].ndim
        wrkindex = 0
        self.dats = []
        self.const = const
        self.grid = []
        # Race_check
        self.race_check = True

        if alginp.lang == 'OPSC':
            l = 0
        elif alginp.lang == 'F90':
            l = 1

        indi = []
        for i in range(self.ndim):
            indi = indi + ['x%d' % i]
        indi = ','.join(indi)
        varform = Idx('%s' % (indi), (l, Symbol('nblock', integer=True)))
        indi = ','.join(['i0-1'])
        varform1 = Idx('%s' % (indi), Symbol('nblock', integer=True))
        self.varform = varform

        self.block = Idx('blk', (l, Symbol('nblock', integer=True)))
        self.grid = self.grid + [self.block.upper]
        self.nblocks = alginp.nblocks
        self.multiblock = alginp.multiblock
        self.blockdims = []
        self.halos = []
        self.gridhalo = []
        self.ranges = {}
        self.iterrange = 0
        self.kernel_ind = 0

        if alginp.lang == 'OPSC':
            if self.multiblock:
                raise ValueError('Implement Multi Block code implementation for OPSC')
            code_file = open(BUILD_DIR+'/OPSC_nssolver.cpp', 'w')
            kernel_file = open(BUILD_DIR+'/auto_kernel.h', 'w')
            self.stencils = {}
            self.sten_ind = 0
            self.sten_name = 'stencil%dD_%%d' % self.ndim
            self.blkname = 'auto_block_OPSC'
            self.kername = 'auto_kernel_%d'
            for dim in range(self.ndim):
                temp = Idx('i%d' % dim, Symbol('nx%dp[blk]' % dim, integer=True))
                self.blockdims.append(temp)
                self.grid = self.grid + [temp.upper+1]
        elif alginp.lang == 'F90':
            if self.multiblock:
                raise ValueError('Implement Multi Block code implementation for F90')
            code_file = open(BUILD_DIR+'/F_serial.f90', 'w')
            kernel_file = open(BUILD_DIR+'/subroutines.f90', 'w')
            modulefile = open(BUILD_DIR+'/param_mod.f90', 'w')
            self.module = []
            for dim in range(self.ndim):
                temp = Idx('i%d' % dim, (1, Symbol('nx%dp' % dim, integer=True)))
                self.blockdims.append(temp)
                self.grid = self.grid + [temp.upper]
            self.kername = 'subroutine_%d'
        else:
            raise ValueError('Implement indexing for language  %s' % (alginp.lang))
        perform_algorithm_checks(alginp, self)
        variab = []

        for va in totvar:
            new = Indexed('%s' % str(va), varform)
            variab = variab + [Indexed('%s' % str(va), varform)]
            for eq in range(len(self.inpeq)):
                self.inpeq[eq] = self.inpeq[eq].subs(va, new)
        # prepare the formulas
        formvar = flatten(list(eq.conser for eq in forms))
        formvar = formvar + flatten(list(eq.variables for eq in forms))
        formvar = list(set(formvar))
        self.const = self.const + list(set(flatten(list(eq.const for eq in forms))))
        for va in formvar:
            new = Indexed('%s' % str(va), varform)
            for eq in range(len(self.forms)):
                self.forms[eq] = self.forms[eq].subs(va, new)
        form_dict = equations_to_dict(self.forms)
        # pprint(self.inpeq)
        variab = set(variab)
        substis = {}
        evaluated = []
        sorteval = []
        self.conser = []
        for con in conser:
            new = Indexed('%s' % str(con), varform)
            evaluated = evaluated + [new]
            sorteval = sorteval + [new]
            self.conser = self.conser + [new]
            self.dats = self.dats + [new]
        variab = variab.difference(set(evaluated))
        # Finished preparing all the stuff now prepare the equations, now the evaluations of primitive
        # Variables
        form_evals = []

        for va in variab:
            val = form_dict.get(va)
            count = count_input_equations(va, self.inpeq)
            count = alginp.ceval
            if count >= alginp.ceval:
                if val:
                    form_evals.append(Eq(va, val))
                    self.dats = self.dats + [va]
                else:
                    raise ValueError('I dont know the formula for %s ' % va)
            else:
                substis[va] = val
                # raise ValueError('Implement how to do for count > ceval')
        if alginp.time_sch == 'RK':

            rkloop = Idx('nrk', (l, alginp.time_ooa))
            time_loop = Idx('iter', (l, Symbol('niter', integer=True)))
            a1 = Indexed('a1', Idx('nrk', (l, alginp.time_ooa)))
            a2 = Indexed('a2', Idx('nrk', (l, alginp.time_ooa)))
            rk = Indexed('rk', Idx('nrk', (l, alginp.time_ooa)))
            self.const = self.const + [a1, a2, time_loop.upper+1]
            saveeq = []
            conser_old = []

        form_evals, evaluated = sort_evals(form_evals, evaluated)
        ders = flatten(list(eq.atoms(Derivative) for eq in self.inpeq))
        ders = list(set(ders))
        der_evals = []
        # get the coefficients of finite differences
        derform(self, alginp)
        tempder = []
        for no, der in enumerate(ders):
            tempder = tempder + [der]
            if der.args[0] in evaluated:
                val = der.args[0]
            else:
                val = substis.get(der.args[0])
            if alginp.ceval == 0:
                count = 100
            elif alginp.ceval == 1000:
                count = 999
            else:
                count = count_input_equations(der, self.inpeq)
            order = len(der.args) - 1
            if val:
                if order == 1 and count >= alginp.ceval:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            temp = alginp.wrkarr % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(var, self)
                            der_evals.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                    else:
                        if alginp.time_sch == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), varform)
                            if alginp.const_timestep:
                                time = Symbol('dt')
                                self.const = self.const + [time]
                            else:
                                time = Indexed('dt', varform)
                                self.dats = self.dats + [time]
                            rkfor = (con - old)/(rk*time)
                            conser_old.append(old)
                            self.dats = self.dats + [old]
                            substis[der] = rkfor
                            saveeq.append(Eq(old, con))

                    # substitute RK3 routine now
                elif order == 2 and count >= alginp.ceval:
                    if der.args[1] == der.args[2]:
                        var = der.subs(der.args[0], val)
                        temp = alginp.wrkarr % (wrkindex)
                        wrkindex = wrkindex + 1
                        new = Indexed('%s' % str(temp), varform)
                        self.dats = self.dats + [new]
                        derf = apply_der(var, self)
                        der_evals.append(Eq(new, derf))
                        substis[der] = new
                    else:
                        expr = Derivative(der.args[0], der.args[2])
                        if substis.get(expr):
                            var = der.subs(expr, substis.get(expr))
                            temp = alginp.wrkarr % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(var, self)
                            der_evals.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                        else:
                            temp = alginp.wrkarr % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(expr, self)
                            substis[expr] = new
                            der_evals.append(Eq(new, derf))
                            self.dats = self.dats + [new]
                            # substis[expr] = new
                            var = der.subs(expr, new)
                            temp = alginp.wrkarr % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(var, self)
                            der_evals.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                elif order == 1 and count < alginp.ceval:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            derf = apply_der(var, self)
                            substis[der] = derf
                    else:
                        if alginp.time_sch == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), varform)
                            if alginp.const_timestep:
                                time = Symbol('dt')
                                self.const = self.const + [time]
                            else:
                                time = Indexed('dt', varform)
                                self.dats = self.dats + [time]
                            rkfor = (con - old)/(rk*time)
                            conser_old.append(old)
                            self.dats = self.dats + [old]
                            substis[der] = rkfor
                            saveeq.append(Eq(old, con))
                elif order == 2 and count < alginp.ceval:
                    if der.args[1] == der.args[2]:
                        var = der.subs(der.args[0], val)
                        derf = apply_der(var, self)
                        substis[der] = derf
                    else:
                        expr = Derivative(der.args[0], der.args[2])
                        if substis.get(expr):
                            var = der.subs(expr, substis.get(expr))
                            derf = apply_der(var, self)
                            # easy_debug
                            # derf = derf.simplify()
                            # easy_debug
                            substis[der] = derf
                        else:
                            derf = apply_der(expr, self)
                            var = der.subs(expr, derf)
                            derf = apply_der(var, self)
                            # for atom in list(derf.atoms(Indexed)):
                                # var = der.subs(expr,atom)
                                # derivative = apply_der(var,self)
                                # derf = derf.subs(atom,derivative)
                            # easy_debug
                            # derf = derf.simplify()
                            # easy_debug
                            substis[der] = derf
                else:
                    raise ValueError('No know;')
            elif count >= alginp.ceval:
                temp = alginp.wrkarr % (wrkindex)
                wrkindex = wrkindex + 1
                new = Indexed('%s' % str(temp), varform)
                self.dats = self.dats + [new]
                var = der.args[0]
                form_evals.append(Eq(new, var))
                # substis[var] = new
                var = der.subs(var, new)
                temp = alginp.wrkarr % (wrkindex)
                wrkindex = wrkindex + 1
                new = Indexed('%s' % str(temp), varform)
                derf = apply_der(var, self)
                der_evals.append(Eq(new, derf))
                substis[der] = new
                self.dats = self.dats + [new]
            elif count < alginp.ceval:
                derf = apply_der(der, self)
                substis[der] = derf
            else:
                pprint(var)
                raise ValueError('This is wrong')

        evaluations = {}
        evalind = 0
        form_evals, sorteval = sort_evals(form_evals, sorteval)
        evaluations[evalind] = form_evals
        evalind = evalind+1
        for der in der_evals:
            evaluations[evalind] = der
            evalind = evalind+1
        # final rk equations

        f = open(BUILD_DIR+'/alg.tex', 'w')
        metadata = {"title":"Algorithm", "author":"Satya P Jammy", "institution":"University of Southampton"}
        header, footer = latex_article_header(metadata)
        f.write(header)
        f.write('The equations are\n\n')
        write_latex(f, self.inpeq, varform)
        f.write('The save state equations are\n\n')
        write_latex(f, saveeq, varform)
        f.write('The evaluations performed are\n\n')
        write_latex(f, evaluations, varform)
        f.write('The substitutions in the equations are\n\n')
        write_latex(f, substis, varform)

        for key, value in substis.iteritems():
            for eqno in range(len(self.inpeq)):
                self.inpeq[eqno] = self.inpeq[eqno].xreplace({key: value})

        # Get the final Runge-Kutta update equations
        tempdict = equations_to_dict(saveeq)
        savedic = dict(zip(tempdict.values(), tempdict.keys()))
        # for race_check errors
        # 1. Find the residue,
        # 2. Update conserve
        # 3. Update old conservative so that there will be no race_errors
        if self.race_check:
            final_eqs = {}
            resdue_eq = []
            upd_conser = []
            upd_old_cons = []
            race_ind = 0
        else:
            final_eqs = []
        for eqno in range(len(self.inpeq)):
            lh = self.inpeq[eqno].lhs
            rh = self.inpeq[eqno].rhs
            if self.race_check:
                race_ind = race_ind + 1
                temp = 'race_eq%d' % race_ind
                new = Indexed('%s' % str(temp), varform)
                self.dats = self.dats + [new]
                resdue_eq = resdue_eq + [Eq(new, rh)]
                temp = solve(Eq(lh, 0), self.conser[eqno])
                nrdr = fraction(lh)
                eqn = nrdr[1]*new + temp[0]
                eqn1 = eqn.xreplace({rk: a1})
                upd_conser = upd_conser + [Eq(self.conser[eqno], eqn1)]
                old = savedic.get(self.conser[eqno])
                eqn1 = eqn.xreplace({rk: a2})
                upd_old_cons = upd_old_cons + [Eq(old, eqn1)]
            else:
                eqn = Symbol('eqn%d' % eqno)
                final_eqs = final_eqs + [Eq(eqn, rh)]
                nrdr = fraction(lh)
                temp = solve(Eq(lh, 0), self.conser[eqno])
                eqn = nrdr[1]*eqn + temp[0]
                eqn1 = eqn.xreplace({rk: a1})
                final_eqs = final_eqs + [Eq(self.conser[eqno], eqn1)]
                old = savedic.get(self.conser[eqno])
                eqn1 = eqn.xreplace({rk: a2})
                final_eqs = final_eqs + [Eq(old, eqn1)]
        # Apply filtering to the equations
        if self.race_check:
            race_ind = 0
            final_eqs[race_ind] = resdue_eq
            race_ind = race_ind+1
            final_eqs[race_ind] = upd_conser
            race_ind = race_ind+1
            final_eqs[race_ind] = upd_old_cons
            race_ind = race_ind+1

        f.write('The equations after substitutions are\n\n')
        write_latex(f, self.inpeq, varform)

        f.write('The final rk3 update equations are\n\n')
        write_latex(f, final_eqs, varform)

        if any(alginp.expfilter):
            self.filtereqs = expfiltering(alginp, self, savedic)
            f.write('The filter equations are\n\n')
            for eq in self.filtereqs:
                write_latex(f, eq, varform)

        alg_template = {'a00': 'header', 'a01': 'defdec', 'a02': 'grid', 'a03': 'init', 'a04': 'Bcst'}
        alg_template['a05'] = 'timeloop'
        alg_template['a06'] = 'saveeq'
        alg_template['a07'] = 'innerloop'
        alg_template['a08'] = 'solveeq'
        alg_template['a09'] = 'Bcs'
        alg_template['a10'] = 'endinnerloop'
        alg_template['a11'] = 'afterinnerloop'
        alg_template['a12'] = 'endtimeloop'
        alg_template['a13'] = 'aftertimeloop'
        alg_template['a14'] = 'footer'

        final_alg = {}
        self.const = list(set(self.const))
        self.dats = list(set(self.dats))

        # Writing kernels and kernel calls
        # bcs(self,alginp)
        if alginp.lang == 'OPSC':
            calls = []
            kerns = []
            if evaluations:
                call, kern = OPSC_write_kernel(evaluations, self)  # Evaluations
                calls = calls + call
                kerns = kerns + kern

            call, kern = OPSC_write_kernel(final_eqs, self)  # Evaluations
            calls = calls + call
            kerns = kerns + kern
            final_alg[alg_template.get('a08')] = calls

            call, kern = OPSC_write_kernel(saveeq, self)
            kerns = kerns + kern
            final_alg[alg_template.get('a06')] = call
            init_eq = []
            # TODO this is for TWOD need to change this to 3D some how
            init_eq = init_eq + [parse_expr('Eq(x, dx0)')]
            init_eq = init_eq + [parse_expr('Eq(y, dx1)')]
            init_eq = init_eq + [parse_expr('Eq(u, sin(x)* cos(y))')]
            init_eq = init_eq + [parse_expr('Eq(v, -cos(x)* sin(y))')]
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], self.conser[0] * parse_expr('u', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], self.conser[0] * parse_expr('v', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            write_latex(f, init_eq, varform)
            call, kern = OPSC_write_kernel(init_eq, self)
            kerns = kerns + kern
            final_alg[alg_template.get('a03')] = call
        elif alginp.lang == 'F90':
            calls = []
            kerns = []
            if evaluations:
                call, kern = F90_write_kernel(evaluations, self, alginp)  # Evaluations
                calls = calls + call
                kerns = kerns + kern
            call, kern = F90_write_kernel(final_eqs, self, alginp)  # Evaluations
            calls = calls + call
            kerns = kerns + kern

            final_alg[alg_template.get('a08')] = calls

            call, kern = F90_write_kernel(saveeq, self, alginp)
            kerns = kerns + kern
            final_alg[alg_template.get('a06')] = call
            init_eq = []
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], parse_expr('sin(x)* cos(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], parse_expr('-cos(x)* sin(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            write_latex(f, init_eq, varform)
            call, kern = F90_write_kernel(init_eq, self, alginp)
            kerns = kerns + kern
            final_alg[alg_template.get('a03')] = call

        # The algortihm
       # Commented for checking Fortran subroutine write
        # if any(alginp.expfilter):
            # for eq in self.filtereqs:
                # call, kern = OPSC_write_kernel(eq,self)
                # final_alg[alg_template.get('a11')] = final_alg[alg_template.get('a11')] + [call]
                # kerns = kerns + kern
            # tempeq = []
            # for eq in self.conser:
                # tempeq = tempeq + [Eq(eq,savedic.get(eq))]
            # call, kern = OPSC_write_kernel(tempeq,self)
            # final_alg[alg_template.get('a11')] = final_alg[alg_template.get('a11')] + [call]
            # kerns = kerns + kern
        # TODO
        if alginp.lang == 'OPSC':
            bcs(self, alginp)
        elif alginp.lang == 'F90':
            self.bccall = ['']

        final_alg[alg_template.get('a00')] = header_code(self, alginp)
        final_alg[alg_template.get('a01')] = defdec(self, alginp)
        final_alg[alg_template['a04']] = self.bccall

        final_alg[alg_template['a09']] = self.bccall

        line_comment = {}
        line_comment['OPSC'] = '//'
        line_comment['F90'] = '!'
        lend = {}
        lend['OPSC'] = ';'
        lend['F90'] = ''
        temp = '%sWrite the grid here if required' % line_comment[alginp.lang]
        temp = temp + '\n\n\n' + '%sGrid writing ends here' % line_comment[alginp.lang]
        final_alg[alg_template.get('a02')] = [temp]
        final_alg[alg_template.get('a11')] = ['%s after  after inner time loop' % line_comment[alginp.lang]]

        tloop, tenloop = loop([time_loop], alginp)
        final_alg[alg_template.get('a05')] = tloop
        final_alg[alg_template.get('a12')] = tenloop
        tloop, tenloop = loop([rkloop], alginp)

        final_alg[alg_template.get('a07')] = tloop
        final_alg[alg_template.get('a10')] = tenloop

        # write the coservative variables to files (including halos)
        final_alg[alg_template.get('a13')] = after_time(self, alginp)
        final_alg[alg_template.get('a14')] = footer_code(self, alginp)

        final_alg['kernels'] = kerns

        write_final_code(alg_template, final_alg, code_file, kernel_file, alginp.lang)
        if alginp.lang == 'F90':
            write_module(self.module, modulefile)
            modulefile.close()
        code_file.close()
        kernel_file.close()

        f.write(footer)
        f.close()

        return

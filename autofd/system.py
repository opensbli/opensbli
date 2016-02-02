#!/usr/bin/env python
import os
import subprocess
from sympy import *
from sympy.parsing.sympy_parser import *

import logging
LOG = logging.getLogger(__name__)

from .equations import *
from .algorithm import *
from .codegen_utils import *
from .fortran import *
from .opsc import *
from .latex import LatexWriter

BUILD_DIR = os.getcwd()


class System(object):

    def __init__(self):
        return

    def prepare(self, equations, formulas, algorithm):
        expanded_equations = flatten(list(e.expandedeq for e in equations))
        expanded_formulas = flatten(list(f.expandedeq for f in formulas))
        variables = flatten(list(e.variables for e in equations))
        conser = flatten(list(e.conser for e in equations))
        constants = list(set(flatten(list(e.constants for e in equations))))
        all_variables = list(set(conser + variables))

        self.ndim = equations[0].ndim
        work_array_index = 0
        self.constants = []
        self.dats = []
        self.grid = []
        # Race_check
        self.race_check = True

        if algorithm.language == 'OPSC':
            l = 0
        elif algorithm.language == 'F90':
            l = 1

        indices = []
        for i in range(self.ndim):
            indices += ['x%d' % i]
        indices = ','.join(indices)
        variable_indices = Idx('%s' % (indices), (l, Symbol('nblock', integer=True)))

        # Single-block and multi-block
        self.block = Idx('blk', (l, Symbol('nblock', integer=True)))
        self.grid = self.grid + [self.block.upper]
        self.nblocks = algorithm.nblocks
        self.multiblock = algorithm.multiblock
        self.blockdims = []
        self.halos = []
        self.gridhalo = []
        self.ranges = {}
        self.iterrange = 0
        self.kernel_ind = 0

        if algorithm.language == 'OPSC':
            if self.multiblock:
                raise NotImplementedError("Implement Multi-Block code for OPSC")
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
        elif algorithm.language == 'F90':
            if self.multiblock:
                raise NotImplementedError("Implement Multi-Block code for F90")
            code_file = open(BUILD_DIR+'/F_serial.f90', 'w')
            kernel_file = open(BUILD_DIR+'/subroutines.f90', 'w')
            module_file = open(BUILD_DIR+'/param_mod.f90', 'w')
            self.module = []
            for dim in range(self.ndim):
                temp = Idx('i%d' % dim, (1, Symbol('nx%dp' % dim, integer=True)))
                self.blockdims.append(temp)
                self.grid = self.grid + [temp.upper]
            self.kername = 'subroutine_%d'
        else:
            raise NotImplementedError('Implement indexing for language %s' % (algorithm.language))

        # Check that all the system's properties are set up according to the algorithm.
        self.sanity_check(algorithm)

        # Create new indexed variables and substitute these into the equations in place of the (non-indexed) variables.
        indexed_variables = []
        for variable in all_variables:
            new = Indexed('%s' % str(variable), variable_indices)
            indexed_variables = indexed_variables + [Indexed('%s' % str(variable), variable_indices)]
            for e in range(len(expanded_equations)):
                expanded_equations[e] = expanded_equations[e].subs(variable, new)

        # Prepare the formulas by indexing their variables.
        formula_variables = flatten(list(f.conser for f in formulas))
        formula_variables += flatten(list(f.variables for f in formulas))
        formula_variables = list(set(formula_variables))
        self.constants += list(set(flatten(list(f.constants for f in formulas))))
        for variable in formula_variables:
            new = Indexed('%s' % str(variable), variable_indices)
            for f in range(len(expanded_formulas)):
                expanded_formulas[f] = expanded_formulas[f].subs(variable, new)
        expanded_formulas_dict = equations_to_dict(expanded_formulas)

        # pprint(expanded_equations)
        indexed_variables = set(indexed_variables)
        substis = {}
        evaluated = []
        sorteval = []
        self.conser = []
        for con in conser:
            new = Indexed('%s' % str(con), variable_indices)
            evaluated = evaluated + [new]
            sorteval = sorteval + [new]
            self.conser = self.conser + [new]
            self.dats = self.dats + [new]
        indexed_variables = indexed_variables.difference(set(evaluated))

        # Finished preparing all the stuff now prepare the equations, now the evaluations of primitive variables
        formula_evaluations = []

        for v in indexed_variables:
            val = expanded_formulas_dict.get(v)
            count = variable_count(v, expanded_equations)
            count = algorithm.evaluation_count
            if count >= algorithm.evaluation_count:
                if val:
                    formula_evaluations.append(Eq(v, val))
                    self.dats = self.dats + [v]
                else:
                    raise ValueError('I dont know the formula for %s ' % v)
            else:
                substis[v] = val
                # raise ValueError('Implement how to do for count > evaluation_count')

        # Runge-Kutta timestepping scheme
        if algorithm.temporal_scheme == 'RK':
            rkloop = Idx('nrk', (l, algorithm.temporal_order))
            time_loop = Idx('iter', (l, Symbol('niter', integer=True)))
            a1 = Indexed('a1', Idx('nrk', (l, algorithm.temporal_order)))
            a2 = Indexed('a2', Idx('nrk', (l, algorithm.temporal_order)))
            rk = Indexed('rk', Idx('nrk', (l, algorithm.temporal_order)))
            self.constants += [a1, a2, time_loop.upper+1]
            saveeq = []
            conser_old = []

        formula_evaluations, evaluated = sort_evals(formula_evaluations, evaluated)
        derivatives = flatten(list(e.atoms(Derivative) for e in expanded_equations))
        derivatives = list(set(derivatives))
        derivative_evaluations = []

        # Get the coefficients of the finite difference scheme
        derivative_formula(self, algorithm)
        tempder = []
        for no, der in enumerate(derivatives):
            tempder = tempder + [der]
            if der.args[0] in evaluated:
                val = der.args[0]
            else:
                val = substis.get(der.args[0])
            if algorithm.evaluation_count == 0:
                count = 100
            elif algorithm.evaluation_count == 1000:
                count = 999
            else:
                count = variable_count(der, expanded_equations)
            order = len(der.args) - 1
            if val:
                if order == 1 and count >= algorithm.evaluation_count:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(var, self)
                            derivative_evaluations.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                    else:
                        if algorithm.temporal_scheme == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), variable_indices)
                            if algorithm.constant_timestep:
                                time = Symbol('dt')
                                self.constants += [time]
                            else:
                                time = Indexed('dt', variable_indices)
                                self.dats = self.dats + [time]
                            rkfor = (con - old)/(rk*time)
                            conser_old.append(old)
                            self.dats = self.dats + [old]
                            substis[der] = rkfor
                            saveeq.append(Eq(old, con))

                    # substitute RK3 routine now
                elif order == 2 and count >= algorithm.evaluation_count:
                    if der.args[1] == der.args[2]:
                        var = der.subs(der.args[0], val)
                        temp = algorithm.work_array % (work_array_index)
                        work_array_index += 1
                        new = Indexed('%s' % str(temp), variable_indices)
                        self.dats = self.dats + [new]
                        derf = apply_der(var, self)
                        derivative_evaluations.append(Eq(new, derf))
                        substis[der] = new
                    else:
                        expr = Derivative(der.args[0], der.args[2])
                        if substis.get(expr):
                            var = der.subs(expr, substis.get(expr))
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(var, self)
                            derivative_evaluations.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                        else:
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(expr, self)
                            substis[expr] = new
                            derivative_evaluations.append(Eq(new, derf))
                            self.dats = self.dats + [new]
                            # substis[expr] = new
                            var = der.subs(expr, new)
                            temp = algorithm.work_array % (work_array_index)
                            work_array_index += 1
                            new = Indexed('%s' % str(temp), variable_indices)
                            derf = apply_der(var, self)
                            derivative_evaluations.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                elif order == 1 and count < algorithm.evaluation_count:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            derf = apply_der(var, self)
                            substis[der] = derf
                    else:
                        if algorithm.temporal_scheme == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), variable_indices)
                            if algorithm.constant_timestep:
                                time = Symbol('dt')
                                self.constants += [time]
                            else:
                                time = Indexed('dt', variable_indices)
                                self.dats = self.dats + [time]
                            rkfor = (con - old)/(rk*time)
                            conser_old.append(old)
                            self.dats = self.dats + [old]
                            substis[der] = rkfor
                            saveeq.append(Eq(old, con))
                elif order == 2 and count < algorithm.evaluation_count:
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

            elif count >= algorithm.evaluation_count:
                temp = algorithm.work_array % (work_array_index)
                work_array_index += 1
                new = Indexed('%s' % str(temp), variable_indices)
                self.dats = self.dats + [new]
                var = der.args[0]
                formula_evaluations.append(Eq(new, var))
                # substis[var] = new
                var = der.subs(var, new)
                temp = algorithm.work_array % (work_array_index)
                work_array_index += 1
                new = Indexed('%s' % str(temp), variable_indices)
                derf = apply_der(var, self)
                derivative_evaluations.append(Eq(new, derf))
                substis[der] = new
                self.dats = self.dats + [new]
            elif count < algorithm.evaluation_count:
                derf = apply_der(der, self)
                substis[der] = derf
            else:
                pprint(var)
                raise ValueError('This is wrong')

        evaluations = {}
        evaluation_index = 0
        formula_evaluations, sorteval = sort_evals(formula_evaluations, sorteval)
        evaluations[evaluation_index] = formula_evaluations
        evaluation_index += 1
        for der in derivative_evaluations:
            evaluations[evaluation_index] = der
            evaluation_index += 1

        # Write out algorithm in LaTeX form
        latex = LatexWriter()
        latex.open(path=BUILD_DIR+'/algorithm.tex')
        metadata = {"title": "Algorithm", "author": "Satya P Jammy", "institution": "University of Southampton"}
        latex.write_header(metadata)

        latex.write_string("The equations are\n\n")
        latex.write_equations(expanded_equations, variable_indices)
        latex.write_string('The save state equations are\n\n')
        latex.write_equations(saveeq, variable_indices)
        latex.write_string('The evaluations performed are\n\n')
        latex.write_equations(evaluations, variable_indices)
        latex.write_string('The substitutions in the equations are\n\n')
        latex.write_equations(substis, variable_indices)

        for key, value in substis.iteritems():
            for eqno in range(len(expanded_equations)):
                expanded_equations[eqno] = expanded_equations[eqno].xreplace({key: value})

        # Get the final Runge-Kutta update equations
        tempdict = equations_to_dict(saveeq)
        savedic = dict(zip(tempdict.values(), tempdict.keys()))
        # for race_check errors
        # 1. Find the residue,
        # 2. Update conserve
        # 3. Update old conservative so that there will be no race_errors
        if self.race_check:
            final_equations = {}
            resdue_eq = []
            upd_conser = []
            upd_old_cons = []
            race_ind = 0
        else:
            final_equations = []
        for eqno in range(len(expanded_equations)):
            lh = expanded_equations[eqno].lhs
            rh = expanded_equations[eqno].rhs
            if self.race_check:
                race_ind = race_ind + 1
                temp = 'race_eq%d' % race_ind
                new = Indexed('%s' % str(temp), variable_indices)
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
                final_equations = final_equations + [Eq(eqn, rh)]
                nrdr = fraction(lh)
                temp = solve(Eq(lh, 0), self.conser[eqno])
                eqn = nrdr[1]*eqn + temp[0]
                eqn1 = eqn.xreplace({rk: a1})
                final_equations = final_equations + [Eq(self.conser[eqno], eqn1)]
                old = savedic.get(self.conser[eqno])
                eqn1 = eqn.xreplace({rk: a2})
                final_equations = final_equations + [Eq(old, eqn1)]

        # Apply filtering to the equations
        if self.race_check:
            race_ind = 0
            final_equations[race_ind] = resdue_eq
            race_ind = race_ind+1
            final_equations[race_ind] = upd_conser
            race_ind = race_ind+1
            final_equations[race_ind] = upd_old_cons
            race_ind = race_ind+1

        latex.write_string('The equations after substitutions are\n\n')
        latex.write_equations(expanded_equations, variable_indices)

        latex.write_string('The final rk3 update equations are\n\n')
        latex.write_equations(final_equations, variable_indices)

        if any(algorithm.expfilter):
            self.filtereqs = expfiltering(algorithm, self, savedic)
            latex.write_string('The filter equations are\n\n')
            for e in self.filtereqs:
                latex.write_equations(e, variable_indices)

        algorithm_template = {'a00': 'header', 'a01': 'defdec', 'a02': 'grid', 'a03': 'init', 'a04': 'Bcst'}
        algorithm_template['a05'] = 'timeloop'
        algorithm_template['a06'] = 'saveeq'
        algorithm_template['a07'] = 'innerloop'
        algorithm_template['a08'] = 'solveeq'
        algorithm_template['a09'] = 'Bcs'
        algorithm_template['a10'] = 'endinnerloop'
        algorithm_template['a11'] = 'afterinnerloop'
        algorithm_template['a12'] = 'endtimeloop'
        algorithm_template['a13'] = 'aftertimeloop'
        algorithm_template['a14'] = 'footer'

        final_algorithm = {}
        self.constants = list(set(self.constants))
        self.dats = list(set(self.dats))

        # Writing kernels and kernel calls
        # bcs(self,algorithm)
        if algorithm.language == 'OPSC':
            calls = []
            kerns = []
            if evaluations:
                call, kern = OPSC_write_kernel(evaluations, self)  # Evaluations
                calls = calls + call
                kerns = kerns + kern

            call, kern = OPSC_write_kernel(final_equations, self)  # Evaluations
            calls = calls + call
            kerns = kerns + kern
            final_algorithm[algorithm_template.get('a08')] = calls

            call, kern = OPSC_write_kernel(saveeq, self)
            kerns = kerns + kern
            final_algorithm[algorithm_template.get('a06')] = call
            init_eq = []
            # TODO this is for 2D need to change this to 3D some how
            init_eq = init_eq + [parse_expr('Eq(x, dx0)')]
            init_eq = init_eq + [parse_expr('Eq(y, dx1)')]
            init_eq = init_eq + [parse_expr('Eq(u, sin(x)* cos(y))')]
            init_eq = init_eq + [parse_expr('Eq(v, -cos(x)* sin(y))')]
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], self.conser[0] * parse_expr('u', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], self.conser[0] * parse_expr('v', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            latex.write_equations(init_eq, variable_indices)
            call, kern = OPSC_write_kernel(init_eq, self)
            kerns = kerns + kern
            final_algorithm[algorithm_template.get('a03')] = call

        elif algorithm.language == 'F90':
            calls = []
            kerns = []
            if evaluations:
                call, kern = F90_write_kernel(evaluations, self, algorithm)  # Evaluations
                calls = calls + call
                kerns = kerns + kern
            call, kern = F90_write_kernel(final_equations, self, algorithm)  # Evaluations
            calls = calls + call
            kerns = kerns + kern

            final_algorithm[algorithm_template.get('a08')] = calls

            call, kern = F90_write_kernel(saveeq, self, algorithm)
            kerns = kerns + kern
            final_algorithm[algorithm_template.get('a06')] = call
            init_eq = []
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], parse_expr('sin(x)* cos(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], parse_expr('-cos(x)* sin(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            latex.write_equations(init_eq, variable_indices)
            call, kern = F90_write_kernel(init_eq, self, algorithm)
            kerns = kerns + kern
            final_algorithm[algorithm_template.get('a03')] = call

        # The algortihm
       # Commented for checking Fortran subroutine write
        # if any(algorithm.expfilter):
            # for eq in self.filtereqs:
                # call, kern = OPSC_write_kernel(eq,self)
                # final_algorithm[algorithm_template.get('a11')] = final_algorithm[algorithm_template.get('a11')] + [call]
                # kerns = kerns + kern
            # tempeq = []
            # for eq in self.conser:
                # tempeq = tempeq + [Eq(eq,savedic.get(eq))]
            # call, kern = OPSC_write_kernel(tempeq,self)
            # final_algorithm[algorithm_template.get('a11')] = final_algorithm[algorithm_template.get('a11')] + [call]
            # kerns = kerns + kern
        # TODO
        if algorithm.language == 'OPSC':
            bcs(self, algorithm)
        elif algorithm.language == 'F90':
            self.bccall = ['']

        final_algorithm[algorithm_template.get('a00')] = header_code(self, algorithm)
        final_algorithm[algorithm_template.get('a01')] = defdec(self, algorithm)
        final_algorithm[algorithm_template['a04']] = self.bccall
        final_algorithm[algorithm_template['a09']] = self.bccall

        temp = '%s Write the grid here if required' % COMMENT_DELIMITER[algorithm.language]
        temp = temp + '\n\n\n' + '%s Grid writing ends here' % COMMENT_DELIMITER[algorithm.language]
        final_algorithm[algorithm_template.get('a02')] = [temp]
        final_algorithm[algorithm_template.get('a11')] = ['%s after  after inner time loop' % COMMENT_DELIMITER[algorithm.language]]

        tloop, tenloop = loop([time_loop], algorithm)
        final_algorithm[algorithm_template.get('a05')] = tloop
        final_algorithm[algorithm_template.get('a12')] = tenloop
        tloop, tenloop = loop([rkloop], algorithm)

        final_algorithm[algorithm_template.get('a07')] = tloop
        final_algorithm[algorithm_template.get('a10')] = tenloop

        # Write the conservative variables to files (including halos)
        final_algorithm[algorithm_template.get('a13')] = after_time(self, algorithm)
        final_algorithm[algorithm_template.get('a14')] = footer_code(self, algorithm)

        final_algorithm['kernels'] = kerns

        write_final_code(algorithm_template, final_algorithm, code_file, kernel_file, algorithm.language)
        if algorithm.language == 'F90':
            write_module(self.module, module_file)
            module_file.close()
        code_file.close()
        kernel_file.close()

        latex.write_footer()
        latex.close()

        return

    def sanity_check(self, algorithm):
        """ Perform sanity checks on the system's properties. """

        if len(algorithm.bcs) == self.ndim:
            pass
        elif (len(algorithm.bcs) > self.ndim):
            raise ValueError('There are more boundary conditions than the number of dimensions')
        elif (len(algorithm.bcs) < self.ndim):
            raise ValueError('There are less boundary conditions than the number of dimensions')

        if len(algorithm.expfilter) == self.ndim:
            pass
        elif (len(algorithm.expfilter) > self.ndim):
            raise ValueError('There are more options for filter than the number of dimensions')
        elif (len(algorithm.expfilter) < self.ndim):
            raise ValueError('There are less options for filter than the number of dimensions')

        return

    def compile(self, algorithm):
        if algorithm.language == 'OPSC':
            # First translate the generated code using the OPSC translator.
            LOG.debug("Translating OPSC code...")
            exit_code = subprocess.call("python $OPS_INSTALL_PATH/../translator/python/c/ops.py OPSC_nssolver.cpp", shell=True)
            if(exit_code != 0):
                # Something went wrong
                LOG.error("Unable to translate OPSC code. Check OPS_INSTALL_PATH?")
        return

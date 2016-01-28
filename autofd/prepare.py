class PreparedEquations(object):

    def __init__(self, eqs, forms, algorithm):
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

        if algorithm.language == 'OPSC':
            l = 0
        elif algorithm.language == 'F90':
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
            modulefile = open(BUILD_DIR+'/param_mod.f90', 'w')
            self.module = []
            for dim in range(self.ndim):
                temp = Idx('i%d' % dim, (1, Symbol('nx%dp' % dim, integer=True)))
                self.blockdims.append(temp)
                self.grid = self.grid + [temp.upper]
            self.kername = 'subroutine_%d'
        else:
            raise NotImplementedError('Implement indexing for language %s' % (algorithm.language))

        algorithm.sanity_check(self)
        variab = []

        for va in totvar:
            new = Indexed('%s' % str(va), varform)
            variab = variab + [Indexed('%s' % str(va), varform)]
            for eq in range(len(self.inpeq)):
                self.inpeq[eq] = self.inpeq[eq].subs(va, new)

        # Prepare the formulas
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

        # Finished preparing all the stuff now prepare the equations, now the evaluations of primitive variables
        form_evals = []

        for va in variab:
            val = form_dict.get(va)
            count = variable_count(va, self.inpeq)
            count = algorithm.evaluation_count
            if count >= algorithm.evaluation_count:
                if val:
                    form_evals.append(Eq(va, val))
                    self.dats = self.dats + [va]
                else:
                    raise ValueError('I dont know the formula for %s ' % va)
            else:
                substis[va] = val
                # raise ValueError('Implement how to do for count > evaluation_count')

        # Runge-Kutta timestepping scheme
        if algorithm.temporal_scheme == 'RK':
            rkloop = Idx('nrk', (l, algorithm.temporal_order))
            time_loop = Idx('iter', (l, Symbol('niter', integer=True)))
            a1 = Indexed('a1', Idx('nrk', (l, algorithm.temporal_order)))
            a2 = Indexed('a2', Idx('nrk', (l, algorithm.temporal_order)))
            rk = Indexed('rk', Idx('nrk', (l, algorithm.temporal_order)))
            self.const = self.const + [a1, a2, time_loop.upper+1]
            saveeq = []
            conser_old = []

        form_evals, evaluated = sort_evals(form_evals, evaluated)
        ders = flatten(list(eq.atoms(Derivative) for eq in self.inpeq))
        ders = list(set(ders))
        der_evals = []

        # Get the coefficients of the finite difference scheme
        derform(self, algorithm)
        tempder = []
        for no, der in enumerate(ders):
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
                count = variable_count(der, self.inpeq)
            order = len(der.args) - 1
            if val:
                if order == 1 and count >= algorithm.evaluation_count:
                    if der.args[1] != Symbol('t'):
                        if substis.get(der):
                            pass
                        else:
                            var = der.subs(der.args[0], val)
                            temp = algorithm.work_array % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(var, self)
                            der_evals.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                    else:
                        if algorithm.temporal_scheme == 'RK':
                            con = der.args[0]
                            temp = '%s_old' % con.base
                            old = Indexed('%s' % str(temp), varform)
                            if algorithm.constant_timestep:
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
                elif order == 2 and count >= algorithm.evaluation_count:
                    if der.args[1] == der.args[2]:
                        var = der.subs(der.args[0], val)
                        temp = algorithm.work_array % (wrkindex)
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
                            temp = algorithm.work_array % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(var, self)
                            der_evals.append(Eq(new, derf))
                            substis[der] = new
                            self.dats = self.dats + [new]
                        else:
                            temp = algorithm.work_array % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(expr, self)
                            substis[expr] = new
                            der_evals.append(Eq(new, derf))
                            self.dats = self.dats + [new]
                            # substis[expr] = new
                            var = der.subs(expr, new)
                            temp = algorithm.work_array % (wrkindex)
                            wrkindex = wrkindex + 1
                            new = Indexed('%s' % str(temp), varform)
                            derf = apply_der(var, self)
                            der_evals.append(Eq(new, derf))
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
                            old = Indexed('%s' % str(temp), varform)
                            if algorithm.constant_timestep:
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
                temp = algorithm.work_array % (wrkindex)
                wrkindex = wrkindex + 1
                new = Indexed('%s' % str(temp), varform)
                self.dats = self.dats + [new]
                var = der.args[0]
                form_evals.append(Eq(new, var))
                # substis[var] = new
                var = der.subs(var, new)
                temp = algorithm.work_array % (wrkindex)
                wrkindex = wrkindex + 1
                new = Indexed('%s' % str(temp), varform)
                derf = apply_der(var, self)
                der_evals.append(Eq(new, derf))
                substis[der] = new
                self.dats = self.dats + [new]
            elif count < algorithm.evaluation_count:
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
        # Final Runge-Kutta step

        # Write out algorithm in LaTeX form
        latex = LatexWriter()
        latex.open(path=BUILD_DIR+'/algorithm.tex')
        metadata = {"title": "Algorithm", "author": "Satya P Jammy", "institution": "University of Southampton"}
        latex.write_header(metadata)

        latex.write_string("The equations are\n\n")
        latex.write_equations(self.inpeq, varform)
        latex.write_string('The save state equations are\n\n')
        latex.write_equations(saveeq, varform)
        latex.write_string('The evaluations performed are\n\n')
        latex.write_equations(evaluations, varform)
        latex.write_string('The substitutions in the equations are\n\n')
        latex.write_equations(substis, varform)

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

        latex.write_string('The equations after substitutions are\n\n')
        latex.write_equations(self.inpeq, varform)

        latex.write_string('The final rk3 update equations are\n\n')
        latex.write_equations(final_eqs, varform)

        if any(algorithm.expfilter):
            self.filtereqs = expfiltering(algorithm, self, savedic)
            latex.write_string('The filter equations are\n\n')
            for eq in self.filtereqs:
                latex.write_equations(eq, varform)

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
        # bcs(self,algorithm)
        if algorithm.language == 'OPSC':
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
            # TODO this is for 2D need to change this to 3D some how
            init_eq = init_eq + [parse_expr('Eq(x, dx0)')]
            init_eq = init_eq + [parse_expr('Eq(y, dx1)')]
            init_eq = init_eq + [parse_expr('Eq(u, sin(x)* cos(y))')]
            init_eq = init_eq + [parse_expr('Eq(v, -cos(x)* sin(y))')]
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], self.conser[0] * parse_expr('u', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], self.conser[0] * parse_expr('v', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            latex.write_equations(init_eq, varform)
            call, kern = OPSC_write_kernel(init_eq, self)
            kerns = kerns + kern
            final_alg[alg_template.get('a03')] = call

        elif algorithm.language == 'F90':
            calls = []
            kerns = []
            if evaluations:
                call, kern = F90_write_kernel(evaluations, self, algorithm)  # Evaluations
                calls = calls + call
                kerns = kerns + kern
            call, kern = F90_write_kernel(final_eqs, self, algorithm)  # Evaluations
            calls = calls + call
            kerns = kerns + kern

            final_alg[alg_template.get('a08')] = calls

            call, kern = F90_write_kernel(saveeq, self, algorithm)
            kerns = kerns + kern
            final_alg[alg_template.get('a06')] = call
            init_eq = []
            init_eq = init_eq + [Eq(self.conser[0], parse_expr('1.0', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[1], parse_expr('sin(x)* cos(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[2], parse_expr('-cos(x)* sin(y)', evaluate=False))]
            init_eq = init_eq + [Eq(self.conser[3], parse_expr('1 + 0.5*(u**2 + v**2)', evaluate=False))]

            latex.write_equations(init_eq, varform)
            call, kern = F90_write_kernel(init_eq, self, algorithm)
            kerns = kerns + kern
            final_alg[alg_template.get('a03')] = call

        # The algortihm
       # Commented for checking Fortran subroutine write
        # if any(algorithm.expfilter):
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
        if algorithm.language == 'OPSC':
            bcs(self, algorithm)
        elif algorithm.language == 'F90':
            self.bccall = ['']

        final_alg[alg_template.get('a00')] = header_code(self, algorithm)
        final_alg[alg_template.get('a01')] = defdec(self, algorithm)
        final_alg[alg_template['a04']] = self.bccall
        final_alg[alg_template['a09']] = self.bccall

        temp = '%s Write the grid here if required' % COMMENT_DELIMITER[algorithm.language]
        temp = temp + '\n\n\n' + '%s Grid writing ends here' % COMMENT_DELIMITER[algorithm.language]
        final_alg[alg_template.get('a02')] = [temp]
        final_alg[alg_template.get('a11')] = ['%s after  after inner time loop' % COMMENT_DELIMITER[algorithm.language]]

        tloop, tenloop = loop([time_loop], algorithm)
        final_alg[alg_template.get('a05')] = tloop
        final_alg[alg_template.get('a12')] = tenloop
        tloop, tenloop = loop([rkloop], algorithm)

        final_alg[alg_template.get('a07')] = tloop
        final_alg[alg_template.get('a10')] = tenloop

        # Write the conservative variables to files (including halos)
        final_alg[alg_template.get('a13')] = after_time(self, algorithm)
        final_alg[alg_template.get('a14')] = footer_code(self, algorithm)

        final_alg['kernels'] = kerns

        write_final_code(alg_template, final_alg, code_file, kernel_file, algorithm.language)
        if algorithm.language == 'F90':
            write_module(self.module, modulefile)
            modulefile.close()
        code_file.close()
        kernel_file.close()

        latex.write_footer()
        latex.close()

        return

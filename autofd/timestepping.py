from .equations import *
from .kernel import *

class TemporalDiscretisation(object):

    def __init__(self,temporal, grid, const_dt, Spatialsolution):
        if const_dt:
            dt = EinsteinTerm('deltat')
            dt.is_constant = True; dt.is_commutative = True
        else:
            raise NotImplementedError("Varying delta t is not implemented in the code")
        self.nstages = temporal.order
        if temporal.scheme == "Forward" and self.nstages == 1:
            self.coeff = None
            pass
        elif temporal.scheme == "RungeKutta" and self.nstages == 3:
            self.coeff = RungeKutta(self.nstages)
            pass
        else:
            raise ValueError("Only 1st order Forward or RungeKutta 3 order temporal schemes are allowed")
        eqs = []
        # Any computations at the start of the time step generally save equations
        self.start_computations = []
        self.computations = []
        self.conservative = []
        # as of now no end computations. This will be updated in the diagnostics
        self.end_computations = None
        out = []
        for soln in Spatialsolution.residual_arrays:
            out.append(self.time_derivative(soln.keys()[0].args[0], dt,soln[soln.keys()[0]], grid))
            self.conservative.append(soln.keys()[0].args[0].base)
        if self.nstages !=1:
            start = [o[-1] for o in out]
            range_ofevaluation = [tuple([0+grid.halos[i][0],s+grid.halos[i][1]]) for i,s in enumerate(grid.shape)]
            self.start_computations.append(Kernel(start,range_ofevaluation, "Save equations"))
            range_ofevaluation = [tuple([0,s]) for i,s in enumerate(grid.shape)]
            # these are the update equations of the variables at t + k where k is rk loop
            eqs = [o[0] for o in out]
            self.computations.append(Kernel(eqs,range_ofevaluation, "Rk new (subloop) update"))
            eqs = [o[1] for o in out]
            self.computations.append(Kernel(eqs,range_ofevaluation, "RK old update"))
        else:
            self.start_computations = None
            range_ofevaluation = [tuple([0,s]) for i,s in enumerate(grid.shape)]
            self.computations.append(Kernel(out,range_ofevaluation, "Euler update"))

        return
        
    def time_derivative(self, fn, dt, residual, grid):
        out = fn
        if self.nstages == 1:
            eqn = Eq(fn, fn + dt*residual, evaluate=False)
        elif self.nstages == 3:
            old = grid.work_array('%s_old'%fn.base)
            eqn_fn = Eq(fn, old + self.coeff.new*residual, evaluate=False)
            eqn_old = Eq(old, old + self.coeff.old*residual, evaluate=False)
            saveeq = Eq(old, fn)
            eqn = [eqn_fn,eqn_old, saveeq]
        return eqn

class RungeKutta():
    def __init__(self, order):
        self.stage = Symbol('stage', integer=True)
        self.old = IndexedBase('rkold')
        self.old.is_grid = False
        self.old.is_constant = True
        self.old.ranges = order
        self.new = IndexedBase('rknew')
        self.new.is_grid = False
        self.new.is_constant = True
        self.new.ranges = order
        self.old = self.old[self.stage]
        self.new = self.new[self.stage]
        if order == 3:
            self.coeffs = {}
            self.coeffs[self.old] = [-1, 2, -1]
            self.coeffs[self.new] = [-1, 4, -6, 4, -1]
        return

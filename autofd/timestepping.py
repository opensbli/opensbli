from .equations import *
from .kernel import *

class TemporalDiscretisation(object):

    """ Perform a temporal discretisation of the equations on the numerical grid of solution points. """

    def __init__(self, temporal, grid, const_dt, spatial_discretisation):
    
        if const_dt:
            dt = EinsteinTerm('deltat')
            dt.is_constant = True
            dt.is_commutative = True
        else:
            raise NotImplementedError("Varying delta t is not implemented in the code.")
            
        self.nstages = temporal.order
        if temporal.scheme == "Forward" and self.nstages == 1:
            self.coeff = None
            pass
        elif temporal.scheme == "RungeKutta" and self.nstages == 3:
            self.coeff = RungeKutta(self.nstages)
            pass
        else:
            raise ValueError("Only first-order Forward or third-order Runge-Kutta temporal discretisation schemes are allowed.")
        
        # Any computations at the start of the time-step. Generally these are the save equations.
        self.start_computations = []
        self.computations = []
        self.conservative = []
        # As of now, no end computations. This will be updated in the diagnostics.
        self.end_computations = None
        
        out = []
        for residual in spatial_discretisation.residual_arrays:
            out.append(self.time_derivative(residual.keys()[0].args[0], dt, residual[residual.keys()[0]], grid))
            self.conservative.append(residual.keys()[0].args[0].base)
            
        if self.nstages != 1:
            # The 'save' equations
            start = [o[-1] for o in out]
            range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
            self.start_computations.append(Kernel(start, range_of_evaluation, "Save equations"))

            # The 'update' equations of the variables at time 't + k', where k is the Runge-Kutta loop iteration.
            range_of_evaluation = [tuple([0, s]) for i, s in enumerate(grid.shape)]
            equations = [o[0] for o in out]
            self.computations.append(Kernel(equations, range_of_evaluation, "RK new (subloop) update"))
            equations = [o[1] for o in out]
            self.computations.append(Kernel(equations, range_of_evaluation, "RK old update"))
        else:
            self.start_computations = None
            range_of_evaluation = [tuple([0, s]) for i, s in enumerate(grid.shape)]
            self.computations.append(Kernel(out, range_of_evaluation, "Euler update"))

        return
        
    def time_derivative(self, fn, dt, residual, grid):
        """ Return the equations used to advance the model equations forward in time. """

        if self.nstages == 1:
            eqn = Eq(fn, fn + dt*residual, evaluate=False)
        elif self.nstages == 3:
            old = grid.work_array('%s_old' % fn.base)
            eqn_fn = Eq(fn, old + self.coeff.new*residual, evaluate=False)
            eqn_old = Eq(old, old + self.coeff.old*residual, evaluate=False)
            save_equation = Eq(old, fn)
            eqn = [eqn_fn, eqn_old, save_equation]
        return eqn


class RungeKutta(object):
    
    """ Runge-Kutta time-stepping scheme. """

    def __init__(self, order):
        """ Set up the Runge-Kutta stages and the coefficients. """
    
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

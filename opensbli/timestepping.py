from .equations import EinsteinTerm
from .scheme import *
from .kernel import *


class TemporalDiscretisation(object):

    """ Perform a temporal discretisation of the equations on the numerical grid of solution points. """

    def __init__(self, temporal_scheme, grid, const_dt, spatial_discretisation):
    
        # Time-stepping scheme
        self.scheme = temporal_scheme
    
        # Constant or variable time-step
        if const_dt:
            dt = EinsteinTerm('deltat')
            dt.is_constant = True
            dt.is_commutative = True
        else:
            raise NotImplementedError("Varying delta t is not implemented in the code.")
        
        # Coefficients required for the time-stepping scheme, from the relevant Butcher tableau.
        self.nstages = temporal_scheme.order
        if isinstance(temporal_scheme, ForwardEuler) and self.nstages == 1:
            self.coeff = None
        elif isinstance(temporal_scheme, RungeKutta) and self.nstages == 3:
            self.coeff = self.scheme.get_coefficients()
        else:
            raise ValueError("Only first-order Forward or third-order Runge-Kutta temporal discretisation schemes are allowed.")
        
        # Start computations: Any computations at the start of the time-step. Generally these are the 'save' equations.
        self.start_computations = []
        
        # Computations: The main computations to perform during each time-step.
        self.computations = []
        
        # The conservative variables (e.g. rho, rhou0, rhou1, rhoE in the Navier-Stokes equations) that should be updated.
        self.conservative = []
        
        # End computations: As of now, no end computations.
        # TODO: This will be updated in the diagnostics.
        self.end_computations = None
        
        # The residual arrays that contain the change in the RHS of each equation.
        out = []
        for residual in spatial_discretisation.residual_arrays:
            out.append(self.time_derivative(residual.keys()[0].args[0], dt, residual[residual.keys()[0]], grid))
            self.conservative.append(residual.keys()[0].args[0].base)
        
        # Formulate each step of the time-stepping scheme here as a computational Kernel.
        range_of_evaluation = [tuple([0, s]) for i, s in enumerate(grid.shape)] # Grid point index 0 to nx (or ny or nz)
        if self.nstages != 1:
            # The 'save' equations.
            start = [o[-1] for o in out]
            range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
            self.start_computations.append(Kernel(start, range_of_evaluation, "Save equations"))

            # The 'update' equations of the variables at time 't + k', where k is the Runge-Kutta loop iteration.
            equations = [o[0] for o in out]
            self.computations.append(Kernel(equations, range_of_evaluation, "RK new (subloop) update"))
            equations = [o[1] for o in out]
            self.computations.append(Kernel(equations, range_of_evaluation, "RK old update"))
        else:
            self.start_computations = None
            self.computations.append(Kernel(out, range_of_evaluation, "Euler update"))

        return
        
    def time_derivative(self, fn, dt, residual, grid):
        """ Return the equations used to advance the model equations forward in time. """

        if self.nstages == 1:
            eqn = Eq(fn, fn + dt*residual, evaluate=False)
        elif self.nstages == 3:
            old = grid.work_array('%s_old' % fn.base)
            eqn_fn = Eq(fn, old + self.scheme.new*residual, evaluate=False)
            eqn_old = Eq(old, old + self.scheme.old*residual, evaluate=False)
            save_equation = Eq(old, fn)
            eqn = [eqn_fn, eqn_old, save_equation]
        return eqn


class RungeKutta(Scheme):
    
    """ Runge-Kutta time-stepping scheme. """

    def __init__(self, order):
        """ Set up the Runge-Kutta stages and the coefficients.
        
        :arg int order: The order of accuracy of the scheme.
        """
    
        Scheme.__init__(self, "RungeKutta", order)
    
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

        return
     
    def get_coefficients(self):
        """ Return the coefficients of the Runge-Kutta update equations. 
        
        :returns: A dictionary of (update_equation, coefficients) pairs.
        :rtype: dict
        """
        
        if self.order == 3:
            coeffs = {}
            coeffs[self.old] = [-1, 2, -1]
            coeffs[self.new] = [-1, 4, -6, 4, -1]
        return coeffs
        

class ForwardEuler(Scheme):
    
    """ First-order forward/explicit Euler time-stepping scheme. """

    def __init__(self):
        """ Set up the forward Euler scheme. """
    
        Scheme.__init__(self, "ForwardEuler", 1)

        return

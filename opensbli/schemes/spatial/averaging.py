"""@brief
   @author David J lusher
   @contributors Satya Pramod Jammy
   @details
"""

from opensbli.physical_models.ns_physics import NSphysics
from sympy import sqrt, factor, Rational, pprint
from opensbli.core.opensbliobjects import DataSet, DataSetBase, EinsteinTerm
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.grid import GridVariable
from opensbli.utilities.helperfunctions import increment_dataset



class Averaging(object):
    def __init__(self, locations):
        self.locations = locations
        return

    def get_locations(self, name, direction, block):
        d = DataSetBase(name, block.shape, block.blocknumber)
        loc = d.location
        loc1, loc2 = loc[:], loc[:]
        loc1[direction] = loc[direction] + self.locations[0]
        loc2[direction] = loc[direction] + self.locations[1]
        dset1, dset2 = d[loc1], d[loc2]
        return dset1, dset2


class SimpleAverage(Averaging):
    def __init__(self, locations):
        self.locations = locations
        Averaging.__init__(self, locations)
        print("Simple averaging is being used for the characteristic system.")
        return

    def average(self, functions, direction, name_suffix, block):
        """Performs a simple average.

        :arg functions: List of function (Symbols) to apply averaging on.
        :arg locations: Relative index used for averaging (e.g. [0,1] for i and i+1)
        :arg direction: Axis of the dataset on which the location should be applied.
        :arg name_suffix: Name to be appended to the functions. """
        avg_equations = []
        for f in functions:
            if isinstance(f, EinsteinTerm):
                name = f.get_base()
                a, b = self.get_locations(name, direction, block)
                avg_equations += [OpenSBLIEq(GridVariable('%s_%s' % (name_suffix, name)), factor((a+b)/2))]
        # Average metric terms with simple average
        averaged_metrics = []
        for item in functions:
            if isinstance(item, DataSet):
                name = item.base.simplelabel()
                a = item
                b = increment_dataset(item, direction, 1)
                averaged_metrics += [OpenSBLIEq(GridVariable('%s_%s' % (name_suffix, name)), factor((a+b)/2))]
        return avg_equations + averaged_metrics


class RoeAverage(Averaging):
    def __init__(self, locations, physics=None):
        print("Roe averaging is being used for the characteristic system")
        self.locations = locations
        Averaging.__init__(self, locations)
        self.physics = physics
        return

    def average(self, functions, direction, name_suffix, block):
        self.direction = direction
        evaluations = []
        # Averaged density rho_hat = sqrt(rho_L*rho_R)
        if self.physics:
            physics = self.physics
        else:
            physics = NSphysics(block)
        rho_L, rho_R = self.get_locations('rho', direction, block)
        evaluations += [sqrt(rho_L*rho_R)]
        grid_vars = [GridVariable('%s_%s' % (name_suffix, 'rho'))]
        # Store inverse factor 1/(sqrt(rho_R)+sqrt(rho_L))
        grid_vars += [GridVariable('%s_%s' % (name_suffix, 'inv_rho'))]
        evaluations += [(sqrt(rho_R)+sqrt(rho_L))**(-1)]

        # Average velocity components
        for i in range(block.ndim):
            velocity = 'u%d' % i
            vel_L, vel_R = self.get_locations(velocity, direction, block)
            evaluations += [(sqrt(rho_L)*vel_L+sqrt(rho_R)*vel_R)*grid_vars[1]]
            grid_vars += [GridVariable('%s_%s' % (name_suffix, velocity))]
        # Averaged kinetic energy in terms of u0, u1, u2 grid variables
        KE = factor(Rational(1, 2)*sum([u**2 for u in grid_vars[2:]]))

        # Calcualte enthalpy h = rhoE + P/rho
        P_L, P_R = self.get_locations('p', direction, block)
        rhoE_L, rhoE_R = self.get_locations('rhoE', direction, block)
        H_L = (rhoE_L + P_L)/rho_L
        H_R = (rhoE_R + P_R)/rho_R
        roe_enthalpy = (sqrt(rho_L)*H_L+sqrt(rho_R)*H_R)*grid_vars[1]
        # Average speed of sound
        a_base = block.location_dataset('a')
        roe_a = sqrt((physics.specific_heat_ratio()-1)*(roe_enthalpy - KE))
        evaluations += [roe_a]
        grid_vars += [GridVariable('%s_%s' % (name_suffix, a_base.base.label))]

        # Average metric terms with simple average
        averaged_metrics = []
        for item in functions:
            if isinstance(item, DataSet):
                name = item.base.simplelabel()
                a = item
                b = increment_dataset(item, direction, 1)
                averaged_metrics += [OpenSBLIEq(GridVariable('%s_%s' % (name_suffix, name)), factor((a+b)/2))]
        return [OpenSBLIEq(x, y) for (x, y) in zip(grid_vars, evaluations)] + averaged_metrics

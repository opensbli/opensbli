"""@brief This contains the Physics Object for Navier-Stokes equations.
   @authors Satya Pramod Jammy, David J Lusher
   @contributors
   @details
"""
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliobjects import DataObject, ConstantObject, DataSetBase
from sympy import S, Rational, pprint, Eq
from opensbli.utilities.helperfunctions import dot
from opensbli.core.parsing import EinsteinEquation
import itertools


class Physics(object):
    """A base class defining physical model's that are to be used for
    codegeneration"""


class PhysicsVariable(object):
    blocknumber = None
    shape = None

    def __init__(self, var):
        if not isinstance(var, DataObject):
            raise ValueError("")
        self._variable = var
        self._relation = None
        self._conservative_relation = None
        return

    @property
    def relation(self):
        """Conver the dataobjects to datasets"""
        return self._relation

    @relation.setter
    def relation(self, formula):
        self._relation = formula
        return

    @property
    def conservative_relation(self):
        return self._conservative_relation

    @conservative_relation.setter
    def conservative_relation(self, formula):
        self._conservative_relation = formula

    @property
    def variable(self):
        # returns the dataset of the current variable at the location
        dsetbase = self.datasetbase
        return dsetbase[dsetbase.location]

    @property
    def datasetbase(self):
        return DataSetBase(self._variable, self.shape, self.blocknumber)


class NSphysics(Physics):
    eqobject = EinsteinEquation()

    def __init__(self, block, **settings):
        self.ndim = block.ndim
        PhysicsVariable.blocknumber = block.blocknumber
        PhysicsVariable.shape = block.shape
        self.set_names(**settings)
        return

    def set_names(self, **settings):
        self._density = PhysicsVariable(DataObject('rho'))
        self._density.relation = self._density.variable
        # DataObject('rho')
        m = self.eqobject.expand('rhou_i', self.ndim, "x", substitutions=[], constants=[])
        m2 = [PhysicsVariable(m1) for m1 in m]
        self._momentum = m2  # Expand rhou_i to number of dimensions
        self._totalenergy = PhysicsVariable(DataObject('rhoE'))
        self._pressure = PhysicsVariable(DataObject('p'))
        self._temperature = PhysicsVariable(DataObject('T'))
        self._speed_of_sound = PhysicsVariable(DataObject('a'))
        m = self.eqobject.expand('u_i', self.ndim, "x", substitutions=[], constants=[])
        m2 = [PhysicsVariable(m1) for m1 in m]
        self._velocity = m2  # Expand u_i to number of dimensions
        # set the relations for the velocity
        for d in range(self.ndim):
            self._velocity[d].relation = self.momentum()[d]/self.density()

        # set the relations for the momentum
        for d in range(self.ndim):
            self._momentum[d].relation = self.velocity()[d]*self.density()

        self._Reynoldsnumber = ConstantObject("Re")
        self._CpbyCv = ConstantObject("gama")  # Ratio of specific heats
        self._Prandtlnumber = ConstantObject("Pr")
        self._Mach = ConstantObject("Minf")

        self._speed_of_sound.relation = (self.specific_heat_ratio()*self.pressure()/self.density())

        self._pressure.relation = (self.specific_heat_ratio() - S.One)*(self.total_energy() - Rational(1, 2)*self.density()*(dot(self.velocity(), self.velocity())))
        self._pressure.conservative_relation = (self.specific_heat_ratio() - S.One)*(self.total_energy() - Rational(1, 2)*(dot(self.momentum(), self.momentum()))/self.density())

        self._temperature.relation = (self.total_energy() - Rational(1, 2)*self.density()*(dot(self.velocity(), self.velocity()))) *\
            (self.specific_heat_ratio()*(self.specific_heat_ratio() - S.One)*self.mach_number()**2) / self.density()
        self._temperature.conservative_relation = (self.total_energy() - Rational(1, 2)*(dot(self.momentum(), self.momentum()))/self.density()) *\
            (self.specific_heat_ratio()*(self.specific_heat_ratio() - S.One)*self.mach_number()**2) / self.density()
        return

    def specific_heat_ratio(self):
        return self._CpbyCv

    def mach_number(self):
        return self._Mach

    @property
    def viscosity(self):
        return

    @viscosity.setter
    def viscosity(self, constant=True):
        if constant:
            self._viscosity = ConstantObject('mu')
        else:
            self._viscosity = PhysicsVariable(DataObject('mu'))

    def density(self, relation=False):
        if not relation:
            return self._density.variable
        else:
            return self._density.relation

    def momentum(self, relation=False):
        if not relation:
            return [m.variable for m in self._momentum]
        else:
            return [m.relation for m in self._momentum]

    def total_energy(self, relation=False):
        if not relation:
            return self._totalenergy.variable
        else:
            raise NotImplementedError("")

    def pressure(self, relation=False, conservative=False):
        if not relation:
            return self._pressure.variable
        else:
            if conservative:
                return self._pressure.conservative_relation
            else:
                return self._pressure.relation

    def speed_of_sound(self, relation=False):
        if not relation:
            return self._speed_of_sound.variable
        else:
            return self._speed_of_sound.relation

    def velocity(self, relation=False):
        if not relation:
            return [m.variable for m in self._velocity]
        else:
            return [m.relation for m in self._velocity]

    def temperature(self, relation=False, conservative=False):
        if not relation:
            return self._temperature.variable
        else:
            if conservative:
                return self._temperature.conservative_relation
            else:
                return self._temperature.relation

    def viscosity(self, relation=False):
        if not relation:
            return [m.variable for m in self._momentum]
        else:
            raise NotImplementedError("")


class StatsVariable(object):
    def __init__(self, var):
        if not isinstance(var, DataObject):
            raise ValueError("")
        self._variable = var
        self._relation = None
        self._conservative_relation = None
        return

    @property
    def relation(self):
        """Conver the dataobjects to datasets"""
        return self._relation

    @relation.setter
    def relation(self, formula):
        self._relation = formula
        return

    @property
    def conservative_relation(self):
        return self._conservative_relation

    @conservative_relation.setter
    def conservative_relation(self, formula):
        self._conservative_relation = formula

    @property
    def variable(self):
        # returns the dataset of the current variable at the location
        dsetbase = self.datasetbase
        return dsetbase[dsetbase.location()]

    @property
    def datasetbase(self):
        return DataSetBase(self._variable)


class NSPhysics_Stats(NSphysics, StatsVariable):
    def __init__(self, ndim, **settings):
        NSphysics.__init__(self, ndim)
        self.ndim = ndim
        # Create the components required for stats
        self.create_stat_components(**settings)
        self.niter = ConstantObject('niter', integer=True)
        from sympy import Int
        self.niter.datatype = Int()
        # Equations to initialise stat arrays to zeros
        self.init = self.init_stats()
        # During time loop stat collection equations
        self.collection = self.stat_collection()
        # Reynolds or Favre averaging after the simulation
        self.average = self.favre_average()
        # self.average = self.reynolds_average()
        return

    def create_stat_components(self, **settings):
        # rho mean equation
        self._rho_mean = StatsVariable(DataObject('rhomean'))
        self._rho_mean.relation = self._rho_mean.variable + self._density.relation
        # momentum component means
        self._momentum_means = [StatsVariable(DataObject('rhou%dmean' % i)) for i in range(self.ndim)]
        for no, component in enumerate(self._momentum_means):
            component.relation = component.variable + self._momentum[no].variable
        # Create mean Reynolds stresses
        indices = [i for i in range(self.ndim)]
        stress_components = sorted(set(tuple(sorted(t)) for t in set(itertools.product(indices, indices))), key=lambda element: (element[0], element[1]))
        self._reynolds_stress_means = [StatsVariable(DataObject('rhou%d%dmean' % (i, j))) for (i, j) in stress_components]
        for no, indices in enumerate(stress_components):
            self._reynolds_stress_means[no].relation = self._reynolds_stress_means[no].variable + \
                                                        (self._momentum[indices[0]].variable*self._momentum[indices[1]].variable)/self._density.variable
        return

    def init_stats(self):
        equations = [Eq(self._rho_mean.variable, 0.0)]
        equations += [Eq(mean.variable, 0.0) for mean in self._momentum_means]
        equations += [Eq(mean.variable, 0.0) for mean in self._reynolds_stress_means]
        for eqn in equations:
            pprint(eqn)
        return equations

    def stat_collection(self):
        equations = [Eq(self._rho_mean.variable, self._rho_mean.relation)]
        equations += [Eq(mean.variable, mean.relation) for mean in self._momentum_means]
        equations += [Eq(mean.variable, mean.relation) for mean in self._reynolds_stress_means]
        for eqn in equations:
            pprint(eqn)
        return

    def reynolds_average(self):
        equations = [Eq(self._rho_mean.variable, self._rho_mean.variable/self.niter)]
        equations += [Eq(mean.variable, mean.variable/(self.niter)) for mean in self._momentum_means]
        equations += [Eq(mean.variable, mean.variable/(self.niter)) for mean in self._reynolds_stress_means]

        return

    def favre_average(self):
        equations = [Eq(self._rho_mean.variable, self._rho_mean.variable/self.niter)]
        equations += [Eq(GridVariable('rmean'), self._rho_mean.variable)]
        equations += [Eq(mean.variable, mean.variable/(self.niter*GridVariable('rmean'))) for mean in self._momentum_means]
        equations += [Eq(mean.variable, mean.variable/(self.niter*GridVariable('rmean'))) for mean in self._reynolds_stress_means]
        for eqn in equations:
            pprint(eqn)
        return equations

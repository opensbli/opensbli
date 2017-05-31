"""This contains the Physics Object for Navier-Stokes equations
"""

class Physics(object):
    """A base class defining physical model's that are to be used for 
    codegeneration"""
class PhysicsVariable(object):
    def __init__(self, var):
        if not isinstance(var, DataObject):
            raise ValueError("")
        self.variable = var
        self._relation = None
        return
    @property
    def relation(self):
        """Conver the dataobjects to datasets"""
        return
    @relation.setter
    def relation(self, formula):
        self._relation = formula
        return
    def dataset(self):
        # returns the dataset of the current variable at the location
        return 
    def datasetbase(self):
        # returns the datasetbase object
        return
class NSphysics(Physics):
    def __init__(self, ndim, **settings):
        self.ndim = ndim
        self.set_names(settings)
        return
    def set_names(self, **settings):
        self._density = PhysicsVariable(DataObject('rho'))
        #DataObject('rho')
        self._momentum = # Expand rhou_i to number of dimensions
        self._totalenergy = PhysicsVariable(DataObject('rhoE'))
        self._pressure = PhysicsVariable(DataObject('p'))
        self._Reynoldsnumber = ConstantObject("Re")
        self._CpbyCv = ConstantObject("gama") # Ratio of specific heats
        self._Prandtlnumber = ConstantObject("Pr")
        self._velocity = # Expand u_i to number of dimensions
        return
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
        if relation=False:
            return self._density.variable
        else:
            raise NotImplementedError("")
    def momentum(self, relation=False):
        return
    def total_energy(self, relation=False):
        return
    def pressure(self, relation=False):
        return
    def velocity(self, relation=False):
        return
    def temperature(self, relation=False):
        return
    def viscosity(self, relation=False):
        return
    

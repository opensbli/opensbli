from opensbli.physical_models.ns_physics import *
from opensbli.core import *
ndim = 2
block = SimulationBlock(ndim, block_number = 0)
settings = {}
settings['viscosity'] = 'constant'
ns = NSphysics(ndim, **settings)
print ns.ndim
pprint(ns.density())
pprint(ns.density(relation=True))
pprint(ns.momentum())
pprint(ns.momentum(relation=True))
pprint(ns.total_energy())
pprint(ns.velocity())
pprint(ns.velocity(relation=True))
pprint(ns.temperature())
pprint(ns.pressure())
pprint(ns.pressure(relation=True))
#t = ns.total_energy
#t.relation = ns.mass
#pprint(t.relation)
#pprint(ns.specific_heat_ratio)


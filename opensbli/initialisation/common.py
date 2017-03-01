"""Place holders for various algorithm locations
"""


class BeforeSimulationStarts(object):
    def __init__(self):
        self.number = 0
        return

class AfterSimulationEnds(object):
    """Place holder for the non simulation equations that are to be solved after the time loop
    example, Output to HDF5, any diagnostics
    """
    def __init__(self):
        self.number = 0
        return
class InTheSimulation(object):
    """Place holder for the non simulation equations that are to be solved with in the time loop
    example, Output to HDF5, any diagnostics
    """
    def __init__(self, frequency):
        self.frequency = frequency
        return
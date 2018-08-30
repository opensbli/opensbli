
#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (c) see License file

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

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

    def __init__(self, frequency=False):
        self.frequency = frequency
        return


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


class SimulationDataType(object):
    @staticmethod
    def set_datatype(types):
        SimulationDataType.dtype = types

    @staticmethod
    def dtype():
        return SimulationDataType.dtype

    @staticmethod
    def opsc():
        return SimulationDataType.dtype.opsc()


class DataType(object):
    pass


class Double(DataType):

    @staticmethod
    def opsc():
        return "double"


class FloatC(DataType):
    @staticmethod
    def opsc():
        return "float"


class UserDefined(DataType):
    """ User defined datatype this is either float or double depending on input"""

    def __init__(self):
        return


class Int(DataType):
    @staticmethod
    def opsc():
        return "int"

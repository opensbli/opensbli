
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

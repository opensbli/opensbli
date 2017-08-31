
class SimulationDataType(object):
    @staticmethod
    def set_datatype(types):
        SimulationDataType.dtype = types

    @staticmethod
    def opsc():
        return SimulationDataType.dtype.opsc()


class DataType(object):
    pass


class Double(DataType):
    @staticmethod
    def opsc():
        return "double"


class Float(DataType):
    def __init__(self):
        return


class UserDefined(DataType):
    """ User defined datatype this is either float or double depending on input"""

    def __init__(self):
        return


class Int(DataType):
    @staticmethod
    def opsc():
        return "int"

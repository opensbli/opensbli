#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

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
from opensbli.core.opensbliequations import NonSimulationEquations, Discretisation, Solution, OpenSBLIEquation
from .common import *
from sympy import pprint, flatten

class opensbliIO(object):
    group_number = 0
    @staticmethod
    def increase_io_group_number():
        opensbliIO.group_number += 1
        return

class iohdf5(opensbliIO):
    def  __new__(cls, arrays = None, save_every = None, **kwargs):
        ret = super(iohdf5,cls).__new__(cls)
        ret.order = 0
        ret.group_number = cls.group_number
        cls.increase_io_group_number()
        if arrays:
            ret.add_arrays(arrays)
        else:
            ret.arrays = None
        if kwargs:
            ret.kwargs = {}
            for key in kwargs:
                ret.kwargs[key.lower()] =  kwargs[key].lower()
        else:
            # Default IO type is write to hdf5
            ret.kwargs = {'iotype':"write"}
        ret.algorithm_place = []
        ret.get_alogorithm_location(save_every=save_every)
        return ret

    def get_alogorithm_location(cls, save_every):
        if save_every:
            cls.algorithm_place += [InTheSimulation(save_every)]
        if cls.kwargs['iotype'] == "write":
            cls.algorithm_place += [AfterSimulationEnds()]
        elif cls.kwargs['iotype'] == "read":
            cls.algorithm_place = [BeforeSimulationStarts()]
        else:
            raise ValueError("")
        return

    def add_arrays(cls, arrays):
        cls.arrays = flatten(arrays)
        return

    def write_latex(cls, latex):
        string = ["HDF5 IO type %s on arrays"%(cls.kwargs['iotype'])]
        string += ["%s"%(d) for d in cls.arrays]
        latex.write_string(' '.join(string))
        return
    @property
    def opsc_code(cls):
        code = []
        if cls.kwargs['iotype'] == "write":
            code += cls.hdf5write_opsc_code()
        elif cls.kwargs['iotype'] == "read":
            code += cls.hdf5read_opsc_code()
        else:
            raise ValueError("Cant classify HDF5io")
        return code
    def hdf5write_opsc_code(cls):
        name = "opensbli_output"
        code = ['char name[80];']
        # generate file name
        code += ['sprintf(name, \"opensbli_output_%d.h5\", iter);']
        dataset_write = []
        for ar in cls.arrays:
            block_name = ar.base.blockname
            dataset_write += ['ops_fetch_dat_hdf5_file(%s, name);' % (ar)]

        # generate the block name
        code += ['ops_fetch_block_hdf5_file(%s,name);'%(block_name)] + dataset_write
        return code
    def hdf5read_opsc_code(cls):

        return

class IoGroup():
    pass

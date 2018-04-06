from .common import InTheSimulation, AfterSimulationEnds, BeforeSimulationStarts
from sympy import flatten
from opensbli.utilities.helperfunctions import dataset_attributes, constant_attributes


class opensbliIO(object):
    group_number = 0

    @staticmethod
    def increase_io_group_number():
        opensbliIO.group_number += 1
        return


class iohdf5(opensbliIO):
    def __new__(cls, arrays=None, save_every=None, **kwargs):
        ret = super(iohdf5, cls).__new__(cls)
        ret.order = 0
        ret.group_number = cls.group_number
        cls.increase_io_group_number()
        if kwargs:
            ret.kwargs = {}
            for key in kwargs:
                ret.kwargs[key.lower()] = kwargs[key].lower()
        else:
            # Default IO type is write to hdf5
            ret.kwargs = {'iotype': "write"}
        ret.algorithm_place = []
        ret.get_algorithm_location(save_every=save_every)
        ret.arrays = []
        if arrays:
            ret.add_arrays(arrays)
        return ret

    def get_algorithm_location(cls, save_every):
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
        cls.arrays += flatten(arrays)
        return
    
    def set_read_from_hdf5_arrays(cls, block):
        if cls.kwargs['iotype'] == "read":
            for ar in cls.arrays:
                if str(ar) in block.block_datasets.keys():
                    dset = block.block_datasets[str(ar)]
                    dset.read_from_hdf5 = True
                    block.block_datasets[str(ar)] = dset
                else:
                    block.block_datasets[str(ar)] = ar
                    block.block_datasets[str(ar)].read_from_hdf5 = True
        return

    def write_latex(cls, latex):
        string = ["HDF5 IO type %s on arrays" % (cls.kwargs['iotype'])]
        string += ["%s" % (d) for d in cls.arrays]
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

    def hdf5write_opsc_code(self):
        code = []
        if "name" in self.kwargs:
            if '.h5' in self.kwargs["name"]:
                name = self.kwargs["name"]
            elif '.' in self.kwargs["name"]:
                raise ValueError("")
            else:
                name = self.kwargs["name"] + '.h5'
            if self.dynamic_fname:
                raise ValueError("dynamic fname not allowed ")
            filename = "\"%s\"" % name
        else:
            name = "opensbli_output"
            code += ['char name[80];']
            if self.dynamic_fname:
                code += ['sprintf(name, \"%s_%%06d.h5\", %s);' % (name, self.control_parameter)]
            else:
                code += ['sprintf(name, \"%s.h5\");' % name]
            filename = "name"
        dataset_write = []
        for ar in self.arrays:
            block_name = ar.base.blockname
            dataset_write += ['ops_fetch_dat_hdf5_file(%s, %s);' % (ar, filename)]

        # generate the block name
        code += ['ops_fetch_block_hdf5_file(%s, %s);' % (block_name, filename)] + dataset_write
        return code

    def hdf5read_opsc_code(cls):
        """To keep the abstraction going return nothing for HDF5 OPSC code"""
        return []
    
    @property
    def evaluated_datasets(cls):
        evaluated = set()
        if cls.kwargs['iotype'] == "read":
            evaluated = evaluated.union(set(self.arrays))
        return evaluated


class IoGroup():
    pass

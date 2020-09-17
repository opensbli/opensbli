"""@brief
   @authors Satya Pramod Jammy, David J lusher
   @contributors
   @details
"""

from opensbli.core.opensbliobjects import DataSet, ConstantIndexed, DataObject
import h5py
from opensbli.code_generation.opsc import rc
from sympy import pprint


def get_min_max_halo_values(halos):
    halo_m = []
    halo_p = []
    if not isinstance(halos, ConstantIndexed):
        for direction in range(len(halos)):
            if halos[direction][0]:
                hal = [d.get_halos(0) for d in halos[direction][0]]
                halo_m += [min(hal)]
            else:
                halo_m += [0]
            if halos[direction][1]:
                hal = [d.get_halos(1) for d in halos[direction][1]]
                halo_p += [max(hal)]
            else:
                halo_p += [0]
        return halo_m, halo_p
    else:
        raise ValueError("")


def increment_dataset(expression, direction, value):
    """ Increments an expression containing datasets by the given increment and direction.
    arg: object: expression: A SymPy expression containing datasets to have their indices updated.
    arg: int: direction: The integer direction to apply the increment. (which DataSet axis to apply to)
    arg: int: value: The positive or negative change to apply to the DataSet's index.
    returns: object: expression: The original expression updated to the new DataSet location."""
    d_objs = expression.atoms(DataObject)
    if len(d_objs) > 0:
        raise TypeError("DataObjects found in the increment_dataset function: %s. These should be converted to DataSets." % ', '.join([str(x) for x in d_objs]))
    for dset in expression.atoms(DataSet):
        loc = list(dset.indices)
        loc[direction] = loc[direction] + value
        new_dset = dset.base[loc]
        expression = expression.replace(dset, new_dset)
    return expression


def dot(v1, v2):
    """Performs the dot product of two variables, they can be lists or single values"""
    out = 0
    if isinstance(v1, list):
        if len(v1) == len(v2):
            for i in range(len(v1)):
                out += v1[i]*v2[i]
            return out
        else:
            raise ValueError("")
    else:
        return v1*v2


def get_inverse_deltas(delta):
    """To reduce divisions we create inverse delta, i.e $inv_0 = 1.0/\left(\delta x0 \right)$ and so on..
    """
    # Check if the delta exists in the rationals defined already
    if delta in rc.existing:
        return rc.existing[delta]
    else:
        # Create a new inverse variable
        name = rc.name
        b, exp = delta.as_base_exp()
        rc.name = "inv_%d"
        inv_delta_name = rc.get_next_rational_constant(delta)
        rc.name = name
        return inv_delta_name


def set_hdf5_metadata(dset, halos, npoints, block):
    """ Function to set hdf5 metadata required by OPS to a dataset."""
    if len(halos) != 2:
        raise ValueError("Two halos should be provided for each dimension.")
    for h in halos:
        if len(h) != block.ndim:
            raise ValueError("halos provided for hdf5 output should be of size %d" % block.ndim)
    # The size of negative halos as a list for all dimensions
    d_m = [halos[i][0] for i in range(2)]
    # The size of positive halos as a list for all dimensions
    d_p = [halos[i][1] for i in range(2)]

    dset.attrs.create("d_p", d_p, dtype="int32")
    dset.attrs.create("d_m", d_m, dtype="int32")
    dset.attrs.create("dim", [1], dtype="int32")
    dset.attrs.create("ops_type", u"ops_dat", dtype="S10")
    dset.attrs.create("block_index", [block.blocknumber], dtype="int32")
    dset.attrs.create("base", [0 for i in range(block.ndim)], dtype="int32")
    dset.attrs.create("type", u"double", dtype="S15")
    dset.attrs.create("block", u"%s" % block.blockname, dtype="S25")
    dset.attrs.create("size", npoints, dtype="int32")
    return


def output_hdf5(array, array_name, halos, npoints, block, **kwargs):
    """ Creates an HDF5 file for reading in data to a simulation,
    sets the metadata required by the OPS library. """
    if not isinstance(array, list):
        array = [array]
    if not isinstance(array_name, list):
        array_name = [array_name]
    assert len(array) == len(array_name)
    if 'filename' in kwargs.keys():
        fname = kwargs['filename']
    else:
        fname = "data.h5"
    with h5py.File(fname, 'w') as hf:
        # Create a group
        g1 = hf.create_group(block.blockname)
        # Loop over all the dataset inputs and write to the hdf5 file
        for ar, name in zip(array, array_name):
            g1.attrs.create("dims", [block.ndim], dtype="int32")
            g1.attrs.create("ops_type", u"ops_block", dtype="S9")
            g1.attrs.create("index", [block.blocknumber], dtype="int32")
            block_dset_name = block.location_dataset(name).base
            dset = g1.create_dataset('%s' % (block_dset_name), data=ar)
            set_hdf5_metadata(dset, halos, npoints, block)
    return


def substitute_simulation_parameters(constants, values, simulation_name='opensbli'):
    """ Function to substitute user provided numerical values for constants
    defined in the simulation.

    :arg list constants: List of strings, one for each input constant in the simulation.
    :arg list values: Numerical values corresponding to the strings in the constants list."""
    file_path = "./%s.cpp" % simulation_name
    substitutions = dict(zip(constants, values))
    print("Constant simulation values:")
    pprint(substitutions)
    with open(file_path) as f:
        s = f.read()
    with open(file_path, 'w') as f:
        for const, value in substitutions.items():
            old_str = const + '=Input;'
            if old_str in s:
                new_str = const + ' = %s' % value + ';'
                s = s.replace(old_str, new_str)
        f.write(s)
    return


def dataset_attributes(dset):
    """
    Move to datasetbase? Should we??
    """
    dset.block_number = None
    dset.read_from_hdf5 = False
    dset.dtype = None
    dset.size = None
    dset.halo_ranges = None
    dset.block_name = None
    return dset


def constant_attributes(const):
    const.is_input = True
    const.dtype = None
    const.value = None
    return const


def print_iteration_ops(simulation_name='opensbli', every=250, NaN_check=None):
    """ Prints the iteration number to standard output. If an array name is passed to NaNcheck
    then the OPS NaN_check is also called. Requires OPS versions since 01/03/2019."""
    file_path = "./%s.cpp" % simulation_name
    with open(file_path) as f:
        lines = f.readlines()
    for no, line in enumerate(lines):
        check_string = "int iter=0;"
        if check_string in line:
            lines[no+1] = lines[no+1] + """if(fmod(iter+1, %d) == 0){
        ops_printf("Iteration is %%d\\n", iter+1); """ % every
            if NaN_check is not None:
                lines[no+1] = lines[no+1] + """
        ops_NaNcheck(%s);\n}\n""" % NaN_check
            else:
                lines[no+1] = lines[no+1] + """\n}\n"""
    with open(file_path, 'w') as f:
        f.write(''.join(lines))
    return

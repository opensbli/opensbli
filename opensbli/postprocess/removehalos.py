"""@brief
   @author Satya Pramod Jammy
   @contributors David J Lusher
   @details
"""

from __future__ import division
import h5py
import argparse

# TODO V2: Documentation on this.


def read_dataset(openname, dataset):
    d_m = openname["%s" % (dataset)].attrs['d_m']
    size = openname["%s" % (dataset)].shape
    read_start = [abs(d) for d in d_m]
    read_end = [s-abs(d) for d, s in zip(d_m, size)]
    if len(read_end) == 1:
        read_data = group["%s" % (dataset)].value[read_start[0]:read_end[0]]
    elif len(read_end) == 2:
        read_data = openname["%s" % (dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1]]
    elif len(read_end) == 3:
        read_data = openname["%s" % (dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
    else:
        raise NotImplementedError("")
    return read_data


def strip_halos(fname_to_read, fname_to_write):
    opensbli_file = h5py.File(fname_to_read, 'r')
    block_name1 = opensbli_file.keys()[0]
    group_block = opensbli_file[block_name1]
    output_opensbli = h5py.File(fname_to_write, 'w')
    # TODO create the same structure of the input file
    for key in group_block.keys():
        data_without_halos = read_dataset(group_block, key)
        output_opensbli.create_dataset("%s" % (key), data=data_without_halos)
    output_opensbli.close()
    return


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user, Teo paths should be provided
    parser = argparse.ArgumentParser(prog="pat")
    parser.add_argument("input_path", help="Path of the HDF5 file written out from OpenSBLI inlcuding the file name", action="store", type=str)
    args = parser.parse_args()
    print "Processing HDF5 from the path %s" % args.input_path
    a = args.input_path.split('/')
    if '.h5' not in a[-1]:
        raise ValueError("Provide the HDF5 file with .h5 extension")
    h5name_output = a[-1].split('.')[0]
    fname_to_write = '/'.join(a[:-1]+['%s_pp.h5' % h5name_output])
    print "Output for post processing will be %s" % fname_to_write
    strip_halos(args.input_path, fname_to_write)

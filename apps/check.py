"""Check the OpenSBLI output file has not become unstable."""

import numpy
import h5py
import argparse

def check(fname):
    opensbli_output_file = h5py.File(fname, 'r')
    # Loop over blocks 
    for block in opensbli_output_file.keys():
        group_block =  opensbli_output_file[block]
        # Loop over all the datasets on block
        for key in group_block.keys():
            # Check if the dataset is NaN
            if numpy.isnan(group_block[key].value).any():
                print("Dataset %s on block %s has  NaN" % (key, block))
            else:
                print("Dataset %s on block %s is OK" % (key, block))

if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="pat")
    parser.add_argument("input_path", help="Path of the HDF5 file written out from OpenSBLI inlcuding the file name", action="store", type=str)
    args = parser.parse_args()
    print("Checking HDF5 from the path %s" % args.input_path)
    check(args.input_path)


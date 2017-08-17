# Convert the output from OpenSBLI to the domain only by removing the halo points
# Author Satya P Jammy, 2017 
# Requires numpy and h5py

from __future__ import division
import numpy as np
import h5py
import argparse


def read_dataset(openname, dataset):
    d_m = openname["%s"%(dataset)].attrs['d_m']
    d_p = openname["%s"%(dataset)].attrs['d_p']
    size = openname["%s"%(dataset)].attrs['size']
    read_start = [abs(d) for d in d_m]
    read_end = [abs(d)+s for d,s in zip(d_m,size)]
    if len(read_end) == 2:
        read_data = openname["%s"%(dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1]]
    elif len(read_end) == 3:
        read_data = openname["%s"%(dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
    else:
        raise NotImplementedError("")
    return read_data

def strip_halos(fname, output_name):
    
    opensbli_file = h5py.File(fname, 'r')
    block_name1 = opensbli_file.keys()[0]
    group_block =  opensbli_file[block_name1]
    output_opensbli = h5py.File(output_name, 'w')
    for key in group_block.keys():
        data_without_halos = read_dataset(group_block, key)
        output_opensbli.create_dataset("%s"%(key), data=data_without_halos)
    
    output_opensbli.close()
    return

if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user, Teo paths should be provided
    parser = argparse.ArgumentParser(prog="pat")
    parser.add_argument("input_path", help="Path of the HDF5 file written out from OpenSBLI inlcuding the file name", action="store", type=str)    
    args = parser.parse_args()
    #fname = 'taylor_green_vortex_500.h5'
    print "Processing HDF5 from the path %s"%args.input_path
    a = args.input_path.split('/')
    if '.h5' not in a[-1]:
        raise ValueError("Provide the HDF5 file with .h5 extension")
    h5name_output = a[-1].split('.')[0]
    output_name = '/'.join(a[:-1]+['%s_pp.h5'%h5name_output])
    print "Output for post processing will be %s"%output_name
    #exit()
    strip_halos(args.input_path, output_name)

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
    size = openname["%s"%(dataset)].shape
    read_start = [abs(d) for d in d_m]
    read_end = [s-abs(d) for d,s in zip(d_m,size)]
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

class Smesh(object):
    coord_names = ["X", "Y", "Z"]
    def __init__(self, dims, size, fname, coordinates=True):
        self.ndim = dims
        self.size = size
        self.node_size_str = " ".join([str(s) for s in self.size])
        self.cell_size = [s-1 for s in size]
        self.fname = fname
        if coordinates:
            self.coordinate = ["x%d_B0"%d for d in range(self.ndim)]
        else:
            self.coordinate = coordinates
        return
    def topology(self):
        return "<Topology TopologyType=\"%dDSMesh\" NumberOfElements=\"%s\"/>" %(self.ndim, self.node_size_str)
    @property
    def dataitem_node(self):
        return """<DataItem Dimensions=\"%s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n %s:/%s\n</DataItem>\n""" 
    def attribute_node(self, attribute):
        attr = """<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n"""%(attribute)
        attr += self.dataitem_node%(self.node_size_str, self.fname, attribute) + "</Attribute>\n"
        return attr 
    def coordinate_read(self):
        geom = "<Geometry GeometryType=\"%s\">\n"%("_".join(self.coord_names[0:self.ndim]))
        reading = geom
        for d in range(self.ndim):
            reading += self.dataitem_node%(self.node_size_str, self.fname, self.coordinate[d])
        reading += "</Geometry>\n"
        return reading

def write_xdmf(output_name):
    opensbli_file = h5py.File(output_name, 'r')
    block_name1 = opensbli_file.keys()[0]
    size =  opensbli_file[block_name1].shape # size of the first data set
    mesh = Smesh(len(size), size, output_name)
    file_write = """<?xml version=\"1.0\" ?>
<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>
<Xdmf Version=\"2.0\">
<Domain>
    <Grid Name=\"mesh1\" GridType=\"Uniform\">
     """+ mesh.topology()+ mesh.coordinate_read() + ''.join([mesh.attribute_node(k) for k in opensbli_file.keys()])+"""   </Grid>\n</Domain>\n</Xdmf>\n"""
    print file_write
    return
    
if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user, Teo paths should be provided
    parser = argparse.ArgumentParser(prog="pat")
    parser.add_argument("input_path", help="Path of the HDF5 file written out from OpenSBLI inlcuding the file name", action="store", type=str)
    args = parser.parse_args()
    #print "Processing HDF5 from the path %s"%args.input_path
    a = args.input_path.split('/')
    if '.h5' not in a[-1]:
        raise ValueError("Provide the HDF5 file with .h5 extension")
    h5name_output = a[-1].split('.')[0]
    output_name = '/'.join(a[:-1]+['%s_pp.h5'%h5name_output])
    #print "Output for post processing will be %s"%output_name
    strip_halos(args.input_path, output_name)
    write_xdmf(output_name)

from math import *
import argparse
import os, os.path

def write_xdmf_file(path_to_hdf5_file, number_of_points, number_of_halos):
    # Number of halos at each end
    halos = 2

    # Number of grid points (including halo points at both ends).
    Nx = number_of_points + number_of_halos*2
    Ny = Nx
    Nz = Nx

    hdf5_filename = os.path.basename(path_to_hdf5_file)
    filename = hdf5_filename.split(".")[0] + '.xmf'
    f = open(filename, 'w')

    # XDMF header
    f.write('''<?xml version="1.0" ?>
    <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
    <Xdmf Version="2.0">
    <Domain>
    ''')

    rhou0 = '%s:/opensbliblock00/rhou0_B0' % path_to_hdf5_file
    rhou1 = '%s:/opensbliblock00/rhou1_B0' % path_to_hdf5_file
    rhou2 = '%s:/opensbliblock00/rhou2_B0' % path_to_hdf5_file
    rho = '%s:/opensbliblock00/rho_B0' % path_to_hdf5_file
    rhoE = '%s:/opensbliblock00/rhoE_B0' % path_to_hdf5_file

    origin_x = 0
    origin_y = 0
    origin_z = 0

    dx = 4.0*pi/Nx
    dy = 2.0/Ny
    dz = 2.0*pi/Nz

    # Grid
    f.write('''
    <Grid Name="Box" GridType="Uniform">
    <Topology TopologyType="3DCoRectMesh" Dimensions="%d %d %d"/>
    <Geometry GeometryType="ORIGIN_DXDYDZ">
       <DataItem DataType="Float" Dimensions="3" Format="XML">%f %f %f</DataItem>
       <DataItem DataType="Float" Dimensions="3" Format="XML">%f %f %f</DataItem>
    </Geometry>
    '''%(Nx, Ny, Nz, origin_x, origin_y, origin_z, dx, dy, dz))

    # rho
    f.write('''\n
    <Attribute Name="rho" AttributeType="Scalar" Center="Node">
    <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8"
    Format="HDF">%s
    </DataItem>
    </Attribute>
    '''%(Nx, Ny, Nz, rho))

    # rhou0
    f.write('''\n
    <Attribute Name="rhou0" AttributeType="Scalar" Center="Node">
    <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8"
    Format="HDF">%s
    </DataItem>
    </Attribute>
    '''%(Nx, Ny, Nz, rhou0))

    # rhou1
    f.write('''\n
    <Attribute Name="rhou1" AttributeType="Scalar" Center="Node">
    <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8"
    Format="HDF">%s
    </DataItem>
    </Attribute>
    '''%(Nx, Ny, Nz, rhou1))


    # rhou2
    f.write('''\n
    <Attribute Name="rhou2" AttributeType="Scalar" Center="Node">
    <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8"
    Format="HDF">%s
    </DataItem>
    </Attribute>
    '''%(Nx, Ny, Nz, rhou2))

    # rhoE
    f.write('''\n
    <Attribute Name="rhoE" AttributeType="Scalar" Center="Node">
    <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8"
    Format="HDF">%s
    </DataItem>
    </Attribute>
    '''%(Nx, Ny, Nz, rhoE))

    # XDMF footer
    f.write('''
    </Grid>
    </Domain>
    </Xdmf>
    ''')
    
    return
    
if(__name__ == "__main__"):
    # Parse the command line arguments
    parser = argparse.ArgumentParser(prog="write_xdmf_file")
    parser.add_argument("path_to_hdf5_file", help="The path to the .h5 file.", action="store", type=str)
    parser.add_argument("number_of_points", help="The number of points, excluding halos (assumed to be the same in all directions).", action="store", type=int)
    parser.add_argument("number_of_halos", help="The number of halo points at each end (assumed to be the same in all directions).", action="store", type=int)
    args = parser.parse_args()

    write_xdmf_file(args.path_to_hdf5_file, args.number_of_points, args.number_of_halos)

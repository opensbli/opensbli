""" Reads in the 2D cylinder grid, adds halos and writes to an HDF5 file.
The HDF5 file should then be given the metadata required by OPS in the cylinder script."""

import h5py
import matplotlib
import numpy as np
from opensbli import *
import os

# Read the grid data
data = np.genfromtxt('./Cyl-grid_Inc.dat')
# nx, ny = data.shape[0,0], data.shape[0,1]
nx, ny = 357, 179
nhalo = 5 # 5 halos on each side of the domain
print("Input grid size: (%d, %d)" % (nx, ny))
# Reshape the data into 2D arrays
output_x = data[:,0].reshape(nx,ny)
output_y = data[:,2].reshape(nx,ny)
# Pad the data with zeros in the halos
full_x, full_y = np.zeros((nx+2*nhalo, ny+2*nhalo)), np.zeros((nx+2*nhalo, ny+2*nhalo))
full_x[nhalo:-(nhalo), nhalo:-(nhalo)] = output_x
full_y[nhalo:-(nhalo), nhalo:-(nhalo)] = output_y
# # Copy coordinate data into the halos manually
# Side 0
# Fill out second index in both periodic condition
print("Grid with halos: (%d, %d)" % (full_x.shape[0], full_x.shape[1]))
# X coordinate array
# Left halos
dx1 = np.abs(full_x[6,:] - full_x[5,:])
full_x[4,:] = full_x[5,:] + dx1
full_x[3,:] = full_x[5,:] + 2*dx1
full_x[2,:] = full_x[5,:] + 3*dx1
full_x[1,:] = full_x[5,:] + 4*dx1
full_x[0,:] = full_x[5,:] + 5*dx1
# Right halos
dx2 = np.abs(full_x[nx+4,:] - full_x[nx+3,:])
full_x[nx+5,:] = full_x[nx+3,:] + dx2
full_x[nx+6,:] = full_x[nx+3,:] + 2*dx2
full_x[nx+7,:] = full_x[nx+3,:] + 3*dx2 
full_x[nx+8,:] = full_x[nx+3,:] + 4*dx2 
full_x[nx+9,:] = full_x[nx+3,:] + 5*dx2 
# Bottom halos
full_x[:,4] = full_x[:, ny-1]
full_x[:,3] = full_x[:, ny-2]
full_x[:,2] = full_x[:, ny-3]
full_x[:,1] = full_x[:, ny-4]
full_x[:,0] = full_x[:, ny-5]
# Top halos
full_x[:,ny+5] = full_x[:,5]
full_x[:,ny+6] = full_x[:,6]
full_x[:,ny+7] = full_x[:,7]
full_x[:,ny+8] = full_x[:,8]
full_x[:,ny+9] = full_x[:,9]

# Y coordinate array
# Left halos
# Add constant spacing dy in each of the halos
dy1 = np.abs(full_y[6,:] - full_y[5,:])
full_y[4,:] = -dy1
full_y[3,:] = -2*dy1
full_y[2,:] = -3*dy1
full_y[1,:] = -4*dy1
full_y[0,:] = -5*dy1
# Right halos
dy2 = np.abs(full_y[nx+4,:] - full_y[nx+3,:])
full_y[nx+5,:] = dy2
full_y[nx+6,:] = 2*dy2
full_y[nx+7,:] = 3*dy2 
full_y[nx+8,:] = 4*dy2 
full_y[nx+9,:] = 5*dy2 
# Bottom halos
# DJL 11/04: Haven't changed these yet, they shouldn't matter (in y direction)
full_y[:,4] = full_y[:, ny-1]
full_y[:,3] = full_y[:, ny-2]
full_y[:,2] = full_y[:, ny-3]
full_y[:,1] = full_y[:, ny-4]
full_y[:,0] = full_y[:, ny-5]
# Top halos
full_y[:,ny+5] = full_y[:,5]
full_y[:,ny+6] = full_y[:,6]
full_y[:,ny+7] = full_y[:,7]
full_y[:,ny+8] = full_y[:,8]
full_y[:,ny+9] = full_y[:,9]


#Transpose the data arrangement
full_x, full_y = full_x.T, full_y.T

# Create an output file and write the data
gf = h5py.File('grid.h5', 'w')
gf.create_dataset('x0', data=full_x)
gf.create_dataset('x1', data=full_y)
gf.close()
# Add HDF5 meta-data to the grid file
b = SimulationBlock(2, block_number=0)
print(b.blockname)
print(b.__dict__)
f = h5py.File('grid.h5', 'r')
x0, x1 = f['x0'], f['x1']
npoints = [nx, ny]
halos = [(-nhalo, nhalo), (-nhalo, nhalo)]
arrays, array_names = [x0, x1], ['x0', 'x1']
output_hdf5(arrays, array_names, halos, npoints, b)
f.close()
os.remove('grid.h5')
print("The cylinder grid has been written to the data.h5 file.")

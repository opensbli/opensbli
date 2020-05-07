""" Reads in the 2D cylinder grid, adds halos and writes to an HDF5 file.
The HDF5 file should then be given the metadata required by OPS in the cylinder script."""

import h5py
import matplotlib
import numpy as np

# Read the grid data
data = np.genfromtxt('./Cyl-grid_Inc.dat')
# nx, ny = data.shape[0,0], data.shape[0,1]
nx, ny = 357, 179
nhalo = 5 # 5 halos on each side of the domain
print(nx, ny)
# Reshape the data into 2D arrays
output_x = data[:,0].reshape(nx,ny)
output_y = data[:,2].reshape(nx,ny)
print(output_y.shape)
# Pad the data with zeros in the halos
full_x, full_y = np.zeros((nx+2*nhalo, ny+2*nhalo)), np.zeros((nx+2*nhalo, ny+2*nhalo))
full_x[nhalo:-(nhalo), nhalo:-(nhalo)] = output_x
full_y[nhalo:-(nhalo), nhalo:-(nhalo)] = output_y
# # Copy coordinate data into the halos manually
# Side 0
dx = np.abs((full_x[5,6] - full_x[5,5]))
print(dx)
left_value = full_x[5,5]
right_value = full_x[5,-6]

print(left_value, right_value)

# Fill out second index in both periodic condition
print(full_x.shape)
# Left halos
full_x[4,:] = full_x[nx-1, :]
full_x[3,:] = full_x[nx-2, :]
full_x[2,:] = full_x[nx-3, :]
full_x[1,:] = full_x[nx-4, :]
full_x[0,:] = full_x[nx-5, :]
# Right halos
full_x[nx+0,:] = full_x[5,:]
full_x[nx+1,:] = full_x[6,:]
full_x[nx+2,:] = full_x[7,:]
full_x[nx+3,:] = full_x[8,:]
full_x[nx+4,:] = full_x[9,:]


full_y[4,:] = full_y[nx-1, :]
full_y[3,:] = full_y[nx-2, :]
full_y[2,:] = full_y[nx-3, :]
full_y[1,:] = full_y[nx-4, :]
full_y[0,:] = full_y[nx-5, :]
# Right halos
full_y[nx+0,:] = full_y[5,:]
full_y[nx+1,:] = full_y[6,:]
full_y[nx+2,:] = full_y[7,:]
full_y[nx+3,:] = full_y[8,:]
full_y[nx+4,:] = full_y[9,:]


#Transpose the data arrangement
full_x, full_y = full_x.T, full_y.T

# Create an output file and write the data
gf = h5py.File('grid.h5', 'w')
gf.create_dataset('x0', data=full_x)
gf.create_dataset('x1', data=full_y)
gf.close()

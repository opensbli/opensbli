""" Grid interpolation of restart data from an OpenSBLI output file to a new mesh.
D. Lusher 19/02/2020 with Dr. G. Coleman (NASA Langley).

Usage: python regrid Nx Ny Nz, for the desired number of (Nx, Ny, Nz) grid points for the output."""
import numpy as np
import h5py
import scipy
from scipy import interpolate
import sys
import matplotlib.pyplot as plt
plt.style.use('classic')
import time


class Utilities(object):
    def __init__(self):
        return

    def read_file(self, fname):
        f = h5py.File(fname, 'r')
        group = f["opensbliblock00"]
        return f, group

    def read_dataset(self, group, dataset):
        d_m = group["%s" % (dataset)].attrs['d_m']
        size = group["%s" % (dataset)].shape
        read_start = [abs(d) for d in d_m]
        read_end = [s-abs(d) for d, s in zip(d_m, size)]
        if len(read_end) == 1:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0]]
        elif len(read_end) == 2:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1]]
        elif len(read_end) == 3:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
        else:
            raise NotImplementedError("")
        return read_data

    def extract_data(self, fname, variables):
        """ Extracts the data from the HDF5 file. Data is stored in the self.data dictionary."""
        f, group1 = self.read_file(fname)
        self.data = {}
        for v in variables:
            temp = v+'_B0'
            self.data[v] = self.read_dataset(group1, temp)
        f.close()
        return

    def set_hdf5_metadata(self, dset):
        """ Function to set hdf5 metadata required by OPS to a dataset."""
        # The size of negative halos as a list for all dimensions
        dset.attrs.create("d_p", [self.n_halos for _ in range(self.ndim)], dtype="int32")
        dset.attrs.create("d_m", [-1*self.n_halos for _ in range(self.ndim)], dtype="int32")
        dset.attrs.create("dim", [1], dtype="int32")
        dset.attrs.create("ops_type", u"ops_dat", dtype="S10")
        dset.attrs.create("block_index", [0], dtype="int32")
        dset.attrs.create("base", [0 for i in range(self.ndim)], dtype="int32")
        dset.attrs.create("type", u"double", dtype="S15")
        dset.attrs.create("block", u"%s" % 'opensbliblock00', dtype="S25")
        dset.attrs.create("size", self.n_points, dtype="int32")
        return

    def output_hdf5(self, variables):
        """ Creates an HDF5 file for reading in data to a simulation,
        sets the metadata required by the OPS library. """
        with h5py.File(output_file_name, 'w') as hf:
            # Create a group
            g1 = hf.create_group('opensbliblock00')
            # Loop over all the dataset inputs and write to the hdf5 file
            for var in variables:
                g1.attrs.create("dims", [self.ndim], dtype="int32")
                g1.attrs.create("ops_type", u"ops_block", dtype="S9")
                g1.attrs.create("index", [0], dtype="int32")
                block_dset_name = '%s_B0' % var
                dset = g1.create_dataset('%s' % (block_dset_name), data=self.new_data[var])
                self.set_hdf5_metadata(dset)
        hf.close()
        return


class ReGrid(Utilities):
    def __init__(self, new_grid):
        """ Example usage from the command line:
        python regrid.py Nx Ny Nz, where Nx Ny Nz are integers for the output mesh."""
        self.nx_new, self.ny_new, self.nz_new = new_grid
        # Output grid distribution in (Nx,Ny,Nz)
        self.n_points = [self.nx_new, self.ny_new, self.nz_new]
        self.ndim = 3
        self.Lx, self.Ly, self.Lz = 4*np.pi, 2.0, 4*np.pi/3.0
        self.stretch_factor = 1.7
        self.n_halos = 5
        self.main()
        return

    def generate_new_y_values(self):
        """ Generates the y coordinate values using the hyperbolic tan stretching function."""
        j = np.array([i for i in range(self.ny_new)])
        y = 0.5*self.Ly*(1.0-((np.tanh(self.stretch_factor*(1.0-2.0*(j/(self.ny_new-1.0)))))/(np.tanh(self.stretch_factor))))-1.0
        return y

    def add_halos(self, variables):
        """ Pads the interpolated 3D data arrays with halo points as required by OpenSBLI."""
        nh = self.n_halos
        for key, dset in self.new_data.items():
            new_dset = np.zeros((2*nh+self.nz_new, 2*nh+self.ny_new, 2*nh+self.nx_new))
            new_dset[nh:-nh, nh:-nh, nh:-nh] = dset
            self.new_data[key] = new_dset
        return

    def generate_3D_coordinates(self, x_new, y_new, z_new):
        # Create 3D arrays to hold the (x,y,z) coordinates
        nx, ny, nz = self.nx_new, self.ny_new, self.nz_new
        x_out = np.zeros((nz, ny, nx))
        y_out, z_out = np.zeros_like(x_out), np.zeros_like(x_out)
        # Populate the 3D arrays from the 1D coordinate arrays
        for j in range(ny):
            for k in range(nz):
                x_out[k,j,:] = x_new
        for i in range(nx):
            for k in range(nz):
                y_out[k,:,i] = y_new
        for i in range(nx):
            for j in range(ny):
                z_out[:,j,i] = z_new
        # Store the new 3D coordinate arrays
        self.new_data['x0'] = x_out
        self.new_data['x1'] = y_out
        self.new_data['x2'] = z_out
        return

    def main(self):
        # Start a timer
        start = time.time()
        variables = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
        f, group1 = self.read_file(input_file)
        x_old, y_old, z_old = self.read_dataset(group1, 'x0_B0'), self.read_dataset(group1, 'x1_B0'), self.read_dataset(group1, 'x2_B0')
        nz_old, ny_old, nx_old = y_old.shape
        # Take a line of the old coordinates
        x_old, y_old, z_old = x_old[0,0,:], y_old[0,:,0], z_old[:,0, 0]
        # Generate new x, z coordinates with linear spacing
        x_new, z_new = np.linspace(x_old[0], x_old[-1], self.nx_new), np.linspace(z_old[0], z_old[-1], self.nz_new)
        # y coordinates for the new mesh
        y_new = self.generate_new_y_values()
        # Dictionary to store the new interpolated data for each variable
        self.new_data = {}
        old_flow_data = self.extract_data(input_file, variables)

        print("Old grid (Nx, Ny, Nz): %s" % str((nx_old, ny_old, nz_old)))
        print("New grid (Nx, Ny, Nz): %s" % str((self.nx_new, self.ny_new, self.nz_new)))
        # Loop over each of the variables to interpolate onto the new mesh
        for var_no, var in enumerate(variables):
            print("Interpolating dataset: %s" % var)
            original_data = self.data[var]
            y_extend = np.zeros((nz_old, self.ny_new, nx_old))
            # Extend in y and interpolate over the grid
            print("Performing the interpolation in the y direction.")
            for k in range(nz_old):
                for i in range(nx_old):
                    data_line = original_data[k, :, i]
                    f = interpolate.splrep(y_old, data_line, s=0)
                    interpolated_y_data = interpolate.splev(y_new, f, der=0)
                    y_extend[k, :, i] = interpolated_y_data
            print("Performing the interpolation in the x direction.")
            x_extend = np.zeros((nz_old, self.ny_new, self.nx_new))
            for k in range(nz_old):
                for j in range(self.ny_new):
                    data_line = y_extend[k, j, :]
                    f = interpolate.splrep(x_old, data_line, s=0)
                    interpolated_x_data = interpolate.splev(x_new, f, der=0)
                    x_extend[k, j, :] = interpolated_x_data
            print("Performing the interpolation in the z direction.")
            z_extend = np.zeros((self.nz_new, self.ny_new, self.nx_new))
            for i in range(self.nx_new):
                for j in range(self.ny_new):
                    data_line = x_extend[:, j, i]
                    f = interpolate.splrep(z_old, data_line, s=0)
                    interpolated_z_data = interpolate.splev(z_new, f, der=0)
                    z_extend[:, j, i] = interpolated_z_data

            # Store the final interpolated data for this flow variable
            self.new_data[var] = z_extend

        self.generate_3D_coordinates(x_new, y_new, z_new)
        # Pad the data with zeros for the halos
        self.add_halos(variables)
        # Write the output file to HDF5
        end = time.time()
        print("Total time taken for the interpolations: %.3f seconds." % (end - start))
        start = time.time()
        variables += ['x0', 'x1', 'x2']
        print(variables)
        self.output_hdf5(variables)
        end = time.time()
        print("Total time taken for writing the output HDF5 file: %.3f seconds." % (end - start))
        return


input_file = './opensbli_output.h5'
output_file_name = 'restart.h5'
# Take the input grid points
if len(sys.argv) > 1:
    new_grid = tuple(sys.argv[1:])
    new_grid = tuple([int(i) for i in new_grid])

PS = ReGrid(new_grid)
""" Grid interpolation of restart data from an OpenSBLI output file to a new mesh.
D. Lusher 19/02/2020 with Dr. G. Coleman (NASA Langley). Modified 17/04 to reduce memory usage. Modified 19/04 to not write x, z coordinates to the restart file."""
import numpy as np
import h5py
import scipy
from scipy import interpolate
import sys
import matplotlib.pyplot as plt
plt.style.use('classic')
import time
import gc


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

    def extract_data(self, fname, variable):
        """ Extracts the data from the HDF5 file. Data is stored in the self.data dictionary."""
        f, group1 = self.read_file(fname)
        output = self.read_dataset(group1, variable + '_B0')
        f.close()
        return output

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

    def instantiate_hdf5_file(self):
        """ Creates an HDF5 file for reading in data to a simulation,
        sets the metadata required by the OPS library. """
        with h5py.File(output_file_name, 'w') as hf:
            # Create a group
            g1 = hf.create_group('opensbliblock00')
        hf.close()
        return

    def write_to_hdf5(self, var_name, output_data):
        """ Creates an HDF5 file for reading in data to a simulation,
        sets the metadata required by the OPS library. """
        print("Writing %s_B0 to the HDF5 output file." % var_name)
        with h5py.File(output_file_name, 'a') as hf:
            # Create a group
            g1 = hf['opensbliblock00']
            g1.attrs.create("dims", [self.ndim], dtype="int32")
            g1.attrs.create("ops_type", u"ops_block", dtype="S9")
            g1.attrs.create("index", [0], dtype="int32")
            block_dset_name = '%s_B0' % var_name
            dset = g1.create_dataset('%s' % (block_dset_name), data=output_data)
            self.set_hdf5_metadata(dset)
        del output_data, dset
        gc.collect()
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
        # Create an empty HDF5 file to append the interpolated data to
        self.instantiate_hdf5_file()
        self.main()
        return

    def generate_new_y_values(self):
        """ Generates the y coordinate values using the hyperbolic tan stretching function."""
        j = np.array([i for i in range(self.ny_new)])
        y = 0.5*self.Ly*(1.0-((np.tanh(self.stretch_factor*(1.0-2.0*(j/(self.ny_new-1.0)))))/(np.tanh(self.stretch_factor))))-1.0
        return y

    def add_halos(self, dset):
        """ Pads the interpolated 3D data arrays with halo points as required by OpenSBLI."""
        nh = self.n_halos
        dset_with_halos = np.zeros((2*nh+self.nz_new, 2*nh+self.ny_new, 2*nh+self.nx_new))
        dset_with_halos[nh:-nh, nh:-nh, nh:-nh] = dset
        print("Using %.2gGB of memory to store this dataset with halo points." % (dset_with_halos.nbytes/1.0e9))
        del dset
        gc.collect()
        return dset_with_halos

    def generate_3D_coordinates(self, x_new, y_new, z_new):
        # Create 3D arrays to hold the (x,y,z) coordinates
        nx, ny, nz = self.nx_new, self.ny_new, self.nz_new
        full_dset = np.zeros((nz, ny, nx))
        # Populate the 3D arrays from the 1D coordinate arrays
        # for j in range(ny):
        #     for k in range(nz):
        #         full_dset[k,j,:] = x_new
        # full_dset = self.add_halos(full_dset)
        # self.write_to_hdf5('x0', full_dset)

        full_dset = np.zeros((nz, ny, nx))
        for i in range(nx):
            for k in range(nz):
                full_dset[k, :, i] = y_new
        full_dset = self.add_halos(full_dset)
        self.write_to_hdf5('x1', full_dset)

        # full_dset = np.zeros((nz, ny, nx))
        # for i in range(nx):
        #     for j in range(ny):
        #         full_dset[:,j,i] = z_new
        # full_dset = self.add_halos(full_dset)
        # self.write_to_hdf5('x2', full_dset)

        del full_dset
        gc.collect()
        return

    def create_x_z_coordinates(self, y_old):
        nz_old, ny_old, nx_old = y_old.shape
        dx, dz = self.Lx/nx_old, self.Lz/nz_old
        x, z = np.zeros(nx_old), np.zeros(nz_old)
        for i in range(nx_old):
            x[i] = i*dx
        for k in range(nz_old):
            z[k] = k*dz
        return x, z

    def plot_check(self, original_coord, original, new_coord, new, var):
        # old_location, new_location = int(original_coord.size/2.0), int(new_coord.size/2.0)
        plt.plot(original_coord, original[-1, :, -1], marker='o', markersize=3.5, markevery=5, label='Original data')
        plt.plot(new_coord, new[-1, :, -1], label='Interpolated')
        plt.legend(loc='best')
        plt.savefig('interpolation_result_%s.pdf' % var, bbox_inches='tight')
        plt.clf()
        return

    def main(self):
        # Start a timer
        start = time.time()
        variables = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
        f, group1 = self.read_file(input_file)
        # x_old, y_old, z_old = self.read_dataset(group1, 'x0_B0'), self.read_dataset(group1, 'x1_B0'), self.read_dataset(group1, 'x2_B0')
        y_old = self.read_dataset(group1, 'x1_B0')
        nz_old, ny_old, nx_old = y_old.shape
        x_old, z_old = self.create_x_z_coordinates(y_old)
        # Free some memory
        f.close()
        del group1
        gc.collect()

        # Take a line of the old coordinates
        # x_old, y_old, z_old = x_old[0,0,:], y_old[0,:,0], z_old[:,0, 0]
        y_old = y_old[0, :, 0]
        # Generate new x, z coordinates with linear spacing
        x_new, z_new = np.linspace(x_old[0], x_old[-1], self.nx_new), np.linspace(z_old[0], z_old[-1], self.nz_new)
        # y coordinates for the new mesh
        y_new = self.generate_new_y_values()
        print("Old grid (Nx, Ny, Nz): %s" % str((nx_old, ny_old, nz_old)))
        print("New grid (Nx, Ny, Nz): %s" % str((self.nx_new, self.ny_new, self.nz_new)))

        # Loop over each of the variables to interpolate onto the new mesh
        for var_no, var in enumerate(variables):
            print("Interpolating dataset: %s" % var)
            original_data = self.extract_data(input_file, var)
            y_extend = np.zeros((nz_old, self.ny_new, nx_old))
            # Extend in y and interpolate over the grid
            print("Performing the interpolation in the y direction.")
            for k in range(nz_old):
                for i in range(nx_old):
                    data_line = original_data[k, :, i]
                    f = interpolate.splrep(y_old, data_line, s=0)
                    interpolated_y_data = interpolate.splev(y_new, f, der=0)
                    y_extend[k, :, i] = interpolated_y_data
            del original_data
            gc.collect()
            print("Performing the interpolation in the x direction.")
            x_extend = np.zeros((nz_old, self.ny_new, self.nx_new))
            for k in range(nz_old):
                for j in range(self.ny_new):
                    data_line = y_extend[k, j, :]
                    f = interpolate.splrep(x_old, data_line, s=0)
                    interpolated_x_data = interpolate.splev(x_new, f, der=0)
                    x_extend[k, j, :] = interpolated_x_data
            del y_extend
            gc.collect()
            print("Performing the interpolation in the z direction.")
            z_extend = np.zeros((self.nz_new, self.ny_new, self.nx_new))
            for i in range(self.nx_new):
                for j in range(self.ny_new):
                    data_line = x_extend[:, j, i]
                    f = interpolate.splrep(z_old, data_line, s=0)
                    interpolated_z_data = interpolate.splev(z_new, f, der=0)
                    z_extend[:, j, i] = interpolated_z_data
            # Make some plots
            # self.plot_check(y_old, original_data, y_new, z_extend, var)
            del x_extend
            gc.collect()
            # Add the halos to this variable
            z_extend = self.add_halos(z_extend)
            # Write this variable to the HDF5 file
            self.write_to_hdf5(var, z_extend)
            del z_extend
            gc.collect()

        del x_old, y_old, z_old, interpolated_y_data, interpolated_x_data, interpolated_z_data
        gc.collect()
        # Generate the coordinate arrays and write to the output file
        self.generate_3D_coordinates(x_new, y_new, z_new)
        # Generate 3D arrays for the new (x,y,z) coordinates
        end = time.time()
        print("Total time taken for re-gridding: %.3f minutes." % ((end - start)/60.0))
        return


input_file = './opensbli_output.h5'
output_file_name = 'restart.h5'
# Take the input grid points
if len(sys.argv) > 1:
    new_grid = tuple(sys.argv[1:])
    new_grid = tuple([int(i) for i in new_grid])

PS = ReGrid(new_grid)

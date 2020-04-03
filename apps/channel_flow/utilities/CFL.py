""" CFL calculation of a simulation output from an OpenSBLI HDF5 file.
D. Lusher 19/02/2020 with Dr. G. Coleman (NASA Langley).

Usage: python CFL.py opensbli_output.h5 """
import numpy as np
import h5py
import sys
from scipy import interpolate
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

class Compute_CFL(Utilities):
    def __init__(self, input_fname):
        # Define constants for the simulation
        self.Re, self.Pr, self.gama = 190.71, 0.7, 1.4
        # Taken from M. Carpenter 2N-low storage RK scheme paper
        rk_scheme = 'RK3'
        if rk_scheme is 'RK3':
            print("Approximating the minimum time-step for the 3rd order 3-stage Runge-Kutta scheme.")
            self.CFL = 1.26
        elif rk_scheme is 'RK4':
            print("Approximating the minimum time-step for the 4th order 5-stage Runge-Kutta scheme.")
            self.CFL = 2.43
        else:
            raise ValueError("The explicit Runge-Kutta scheme CFL limit must be set.")
        self.main(input_fname)
        return

    def extract_data(self, fname, variables):
        """ Extracts the data from the HDF5 file. Data is stored in the self.data dictionary."""
        f, group1 = self.read_file(fname)
        # Create a dictionary if it doesn't already exist
        if not hasattr(self, 'data'):
            self.data = {}
        for v in variables:
            temp = v+'_B0'
            self.data[v] = self.read_dataset(group1, temp)
        f.close()
        return

    def perform_average(self, input_variables):
        """ Averages over the two homogeneous directions."""
        n_vars = len(self.data.keys())
        # Write the j indices in the y direction
        output_data = np.array([i+1 for i in range(self.ny)])
        for key in input_variables:
            # if key is 'T_mean' or key is 'u0mean':
            print("Averaging %s" % key)
            # Average over the x direction, exclusing the last plane
            temp_data = self.data[key][:,:,0:-1]
            temp_data = np.mean(temp_data, axis=2)
            temp_data = temp_data[0:-1,:]
            # Average over the z direction, excluding the last plane
            temp_data = np.mean(temp_data, axis=0)
            output_data = np.column_stack([output_data, temp_data])
        return output_data

    def calculate_primitive_variables(self):
        """ Creates the primitive variables from the conservative input."""
        rho = self.data['rho']
        u = self.data['rhou0'] / rho
        v = self.data['rhou1'] / rho
        w = self.data['rhou2'] / rho
        p = (self.gama - 1.0)*(self.data['rhoE'] - 0.5*(u**2 + v**2 + w**2)*rho)
        a = np.sqrt(self.gama*p / rho) # Speed of sound
        return u, v, w, p, a

    def calculate_deltas(self):
        # Get the N-1 grid spacings in x, y, z
        nx, ny, nz = self.nx, self.ny, self.nz
        dx_line = np.diff(self.data['x0'][0,0,:])
        dy_line = np.diff(self.data['x1'][0,:,0])
        dz_line = np.diff(self.data['x2'][:,0,0])
        # Store N grid spacings by fixing the end points and using dx[i] = 0.5*(dx[i+1]+dx[i-1])
        dx_full, dy_full, dz_full = np.zeros(self.nx), np.zeros(self.ny), np.zeros(self.nz)
        # Grid spacings in x
        dx_full[0], dx_full[-1] = dx_line[0], dx_line[-1]
        for i in range(1, self.nx-1):
            dx_full[i] = 0.5*(dx_line[i] + dx_line[i-1])
        # Grid spacings in y 
        dy_full[0], dy_full[-1] = dy_line[0], dy_line[-1]
        for i in range(1, self.ny-1):
            dy_full[i] = 0.5*(dy_line[i] + dy_line[i-1])
        # Grid spacings in z
        dz_full[0], dz_full[-1] = dz_line[0], dz_line[-1]
        for i in range(1, self.nz-1):
            dz_full[i] = 0.5*(dz_line[i] + dz_line[i-1])
        # Copy the dx, dy, dz line profiles to 3D fields to allow array operations
        dx_3D, dy_3D, dz_3D = np.zeros((nz, ny, nx)), np.zeros((nz, ny, nx)), np.zeros((nz, ny, nx))
        for j in range(self.ny):
            for k in range(self.nz):
                dx_3D[k, j, :] = dx_full
        for i in range(self.nx):
            for k in range(self.nz):
                dy_3D[k, :, i] = dy_full
        for i in range(self.nx):
            for j in range(self.ny):
                dz_3D[:, j, i] = dz_full
        return dx_3D, dy_3D, dz_3D

    def compute_CFL(self, rho, u, v, w, p, a, dx, dy, dz):
        # Compute the CFL condition from pg36 of Qinling Li PhD Thesis, Univ Soton.(2003). 
        f1 = np.sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2)
        f2 = 1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2 + 2.0/(dx*dy) + 2.0/(dx*dz) + 2.0/(dy*dz)
        b = np.abs(u)/dx + np.abs(v)/dy + np.abs(w)/dz + f1*a + 4.0*f2*(self.gama/(self.Re*self.Pr))
        return b

    def main(self, input_fname):
        # Variables to extract
        f, group = self.read_file(input_fname)
        input_variables = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
        grid_variables = ['x0', 'x1', 'x2']
        # Read the input file
        start = time.time()
        self.extract_data(input_fname, input_variables + grid_variables)
        end = time.time()
        print("Time taken to read the input file: %.3f seconds.\n" % (end-start))
        start = time.time()
        self.nz, self.ny, self.nx = self.data.values()[0].shape
        print("Performing the CFL calculation on a grid of (Nx,Ny,Nz) = (%d,%d,%d)." % (self.nx, self.ny, self.nz))
        print("Input parameters: Re=%.3f, Pr=%.3f, Cp/Cv=%.3f, CFL=%.3f.\n" % (self.Re, self.Pr, self.gama, self.CFL))
        # Calculate 3D fields of the local grid spacings
        dx_3D, dy_3D, dz_3D = self.calculate_deltas()
        # Calculate primitive quantities
        u, v, w, p, a = self.calculate_primitive_variables()
        # Compute the CFL field in 3D
        rho = self.data['rho']
        b = self.compute_CFL(rho, u, v, w, p, a, dx_3D, dy_3D, dz_3D)

        end = time.time()
        print("Time taken to compute the CFL: %.3f seconds.\n" % (end-start))
        # Plot the y variation of the max time-step
        y_line = self.data['x1'][0,:,0]
        dy_line = dy_3D[0,:,0]
        b_reduce = np.mean(b, axis=2)
        b_reduce = np.mean(b_reduce, axis=0)
        dt = self.CFL/b_reduce
        grid_indices = np.array([j for j in range(self.ny)])
        labels = ["j index", "y coordinate", "b", "max delta t"]
        print("{: >15} {: >15} {: >15} {: >15}".format(*labels))
        for j in grid_indices:
            row = [grid_indices[j], y_line[j], b_reduce[j], dt[j]]
            print("{: >15} {: >15} {: >15} {: >15}".format(*row))
        print("The minimum time-step value (CFL/b_max) over the entire 3D domain is dt=%.6f" % (self.CFL/np.max(b)))
        return

try:
    input_fname = sys.argv[1]
except:
    print("Please provide a filename for the input data as a command line argument.")
    exit()
# Create an instance of the stats class
PS = Compute_CFL(input_fname)

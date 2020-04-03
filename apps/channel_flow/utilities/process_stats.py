""" Spatial averaging of an OpenSBLI HDF5 statistics file. The output is formatted
to be read in by TecPlot. D. Lusher 19/02/2020 with Dr. G. Coleman (NASA Langley).

Usage: python process_stats.py output_file_name.txt case_name."""
import numpy as np
import h5py
import sys

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


class ProcessStats(Utilities):
    def __init__(self, stats_fname):
        # Channel dimensions
        self.Lx, self.Ly, self.Lz = 4*np.pi, 2.0, 4*np.pi/3.0
        self.stretch_factor = 1.7
        self.main(stats_fname)
        return

    def extract_data(self, fname, variables):
        """ Extracts the data from the HDF5 file. Data is stored in the self.data dictionary."""
        f, group1 = self.read_file(fname)
        self.data = {}
        for v in variables:
            temp = v+'_B0'
            self.data[v] = self.read_dataset(group1, temp)
        f.close()
        return

    def generate_grid(self):
        """ Generates the y coordinate values using the hyperbolic tan stretching function."""
        j = np.array([i for i in range(self.ny)])
        y = 0.5*self.Ly*(1.0-((np.tanh(self.stretch_factor*(1.0-2.0*(j/(self.ny-1.0)))))/(np.tanh(self.stretch_factor))))-1.0
        return y

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

    def write_output_data(self, array_names, output_variable_names, output_data):
        """ Write the y profiles to the output dat file."""
        # Write the preamble
        f = open(output_file_name, 'w')
        f.write('VARIABLES = \"%s\"\n' % 'myJ')
        for var in output_variable_names:
            f.write('\"mean_%s\"\n' % var)
        f.write('\"Y\"\n')
        # Write the zone name
        f.write('ZONE T=\"%s\"\n' % zone_name)
        f.write('I=%d\n' % self.ny)
        f.close()
        with open("%s" % output_file_name, "ab") as f:
            np.savetxt(f, output_data, fmt='%1.9E', delimiter=' ')
        return

    def main(self, stats_fname):
        # Variables to extract
        f, group = self.read_file(stats_fname)
        # Variable names in the file
        input_variables = ['u2u2mean', 'u2u1mean', 'u2u0mean', 'u2mean', 'u1u1mean', 'u1u0mean', 'u1mean', 'u0u0mean', 'u0mean', 'rhou2u2mean', 'rhou2u1mean',\
                           'rhou2u0mean', 'rhou2mean', 'rhou1u1mean', 'rhou1u0mean', 'rhou1mean', 'rhou0u0mean', 'rhou0mean', 'rhomean', \
                         'pp_mean', 'p_mean', 'mu_mean', 'a_mean', 'T_mean', 'TT_mean', 'M_mean', 'E_mean', 'vorticity_sq_mean', 'dilatation_sq_mean'] # "rhorho_mean, 'rhoTu1_mean', 'rhoT_mean' Taken out
        # Output names to match the Tecplot example
        variable_names = ["ww", "vw", "uw", "W", "vv", "uv", "V", "uu", "U", "rhoww", "rhovw", "rhouw", "rhow",\
                         "rhovv", "rhouv", "rhov", "rhouu", "rhou", "rho", "pp", "p",\
                         "mu", "a", "T", "TT", "M", "E", "vortsq", "dilatsq"] # "rhorho", "rhoTv", "rhoT" Taken out.
        # Read the statistics file
        self.extract_data(stats_fname, input_variables)
        self.nz, self.ny, self.nx = self.data.values()[0].shape
        # Generate the y coordinates
        y = self.generate_grid()
        # Average over the two periodic dimensions
        output_data = self.perform_average(input_variables)
        # Add the y coordinates to the output data
        output_data = np.column_stack((output_data, y))
        # Write to file
        self.write_output_data(input_variables, variable_names, output_data)
        return

stats_fname = './stats_central_4.h5'
try:
    output_file_name = sys.argv[1]
    zone_name = sys.argv[2]
except:
    print("Please provide a filename for the output as a command line argument.")
    exit()
# Create an instance of the stats class
PS = ProcessStats(stats_fname)

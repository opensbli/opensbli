import fileinput
from shutil import copyfile
""" Script to modify a defdec_data_set.h file to allow a simulation to be restarted. Note the 
initialisation kernel should also be commented out in the main opensbli.cpp file."""

# Arrays to be read in from a restart file, add additional datasets here as required
arrays = ['rho_B0', 'rhou0_B0', 'rhou1_B0', 'rhou2_B0', 'rhoE_B0', 'rho_filt_B0', 'rhou0_filt_B0', 'rhou1_filt_B0', 'rhou2_filt_B0', 'rhoE_filt_B0']
filename = 'defdec_data_set.h'
# Back up the existing declaration file
copyfile('./%s' % filename, './init_defdec.h')

for variable in arrays:
    f = open(filename, 'r')
    file_data = f.read()
    f.close()
    text_to_search = 'ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "%s");' % variable
    replacement_text = 'ops_decl_dat_hdf5(opensbliblock00, 1, "double", "%s", "restart.h5");' % variable
    newdata = file_data.replace(text_to_search, replacement_text)
    f = open(filename, 'w')
    f.write(newdata)
    f.close()

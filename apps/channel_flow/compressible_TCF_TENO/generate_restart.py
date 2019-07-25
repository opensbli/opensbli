import fileinput

arrays = ['rho_B0', 'rhou0_B0', 'rhou1_B0', 'rhou2_B0', 'rhoE_B0', 'x0_B0', 'x1_B0', 'x2_B0']
filename = 'defdec_data_set.h'

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

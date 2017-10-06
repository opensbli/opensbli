from sympy import pprint


def substitute_parameters(simulation_name, constants, values, dsets, hdf5=False):
    file_path = "./%s.cpp" % simulation_name
    substitutions = dict(zip(constants, values))
    pprint(substitutions)
    with open(file_path) as f:
        s = f.read()
    with open(file_path, 'w') as f:
        for const, value in substitutions.iteritems():
            old_str = const + '=Input;'
            if old_str in s:
                new_str = const + ' = %s' % value + ';'
                s = s.replace(old_str, new_str)
        ops_exit = 'ops_exit();'
        if not hdf5:
            if ops_exit in s:
                new_str = ''
                for dset in dsets:
                    new_str += 'ops_print_dat_to_txtfile(%s_B0, "%s.dat");\n' % (dset, dset)
                new_str += ops_exit
        if hdf5:
            if ops_exit in s:
                new_str = 'ops_fetch_block_hdf5_file(%sblock00, "%s.h5");\n' % (simulation_name, simulation_name)
                for dset in dsets:
                    new_str += 'ops_fetch_dat_hdf5_file(%s_B0, "%s.h5");\n' % (dset, simulation_name)
                new_str += ops_exit

        s = s.replace(ops_exit, new_str)
        f.write(s)
    return


if __name__ == "__main__":
    constants = ['gama', 'Minf', 'Pr', 'Re', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2',
                 'Delta0block0', 'Delta1block0', 'Delta2block0', 'SuthT', 'RefT', 'eps', 'TENO_CT', 'Lx1', 'Lx2', 'by', 'bz']
    values = ['1.4', '2.0', '0.72', '950.0', '1.67619431', '1e-2', '0', '600', '500', '150',
              '400.0/(block0np0-1)', '115.0/(block0np1-1)', '50.0/(block0np2-1)', '110.4', '288.0', '1e-15', '1e-5', '115.0', '50.0', '3.0', '3.0']
    simulation_name = 'opensbli'
    dsets = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE', 'x0', 'x1', 'x2']
    substitute_parameters(simulation_name, constants, values, dsets, True)

import os
from sympy import pprint
import numpy as np

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
                if isinstance(value, int):
                    new_str = const + ' = %d;' % value
                elif isinstance(value, float):
                    new_str = const + ' = %f;' % value
                s = s.replace(old_str, new_str)
        # f.write(s)
        # Temporary output dats

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
    constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0', "c0", "c1", "c2"]
    values = [1.0, 1.4, 0.1, 0.7, 0.0025, 5, 64, 64, 64, 4*np.pi/64, 2.0/64, 2*np.pi/64, 0, 0, 0]
    simulation_name = 'opensbli'
    dsets = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
    substitute_parameters(simulation_name, constants, values, dsets, True)

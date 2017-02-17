import os
from sympy import pprint

def substitute_parameters(simulation_name, constants, values, dsets):
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
        if ops_exit in s:
            new_str = ''
            for dset in dsets:
                new_str += 'ops_print_dat_to_txtfile(%s_B0, "%s.dat");\n' % (dset, dset)
            new_str += ops_exit
        s = s.replace(ops_exit, new_str)
        f.write(s)
    return




if __name__ == "__main__":
    constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np_0', 'block0np_1', 'block0np_2', 'Delta_0block_0', 'Delta_1block_0', 'Delta_2block_0']
    values = [1600.0, 1.4, 0.1, 0.71, 3.385*10**-4, 500, 64, 64, 64, 1.0/64, 1.0/64, 1.0/64]
    simulation_name = 'OpenSBLI'
    dsets = ['rho', 'rhou0', 'rhou1', 'rhou2', 'rhoE']
    substitute_parameters(simulation_name, constants, values, dsets)
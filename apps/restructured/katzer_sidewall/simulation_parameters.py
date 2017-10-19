from sympy import pprint

def substitute_parameters(simulation_name, constants, values):
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
        f.write(s)
    return


if __name__ == "__main__":
    constants = ['gama', 'Minf', 'Pr', 'Re', 'Twall', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2',
                 'Delta0block0', 'Delta1block0', 'Delta2block0', 'SuthT', 'RefT', 'eps', 'TENO_CT', 'Lx1', 'Lx2', 'by', 'bz']
    values = ['1.4', '2.0', '0.72', '950.0', '1.67619431', '0.05', '100', '400', '400', '100',
              '400.0/(block0np0-1)', '115.0/(block0np1-1)', '57.5/(block0np2-1)', '110.4', '288.0', '1e-15', '1e-5', '115.0', '57.5', '3.0', '3.0']
    simulation_name = 'opensbli'
    substitute_parameters(simulation_name, constants, values)

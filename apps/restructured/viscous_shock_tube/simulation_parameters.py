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
    constants = ['mu', 'gama', 'Minf', 'Pr', 'Re', 'dt', 'niter', 'block0np0', 'block0np1',
                 'Delta0block0', 'Delta1block0', 'eps', 'TENO_CT']
    values = ['1.0', '1.4', '1.0', '0.73', '200.0', '0.00005', '20000', '600', '300',
              '1.0/(block0np0-1)', '0.5/(block0np1-1)', '1e-15', '1e-7']
    simulation_name = 'opensbli'
    substitute_parameters(simulation_name, constants, values)

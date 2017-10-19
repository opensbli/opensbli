
def substitute_parameters(simulation_name, constants, values):
    file_path = "./%s.cpp" % simulation_name
    substitutions = dict(zip(constants, values))
    print(substitutions)
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
    constants = ['Re', 'gama', 'Minf', 'Pr', 'dt', 'niter', 'block0np0', 'block0np1', 'block0np2', 'Delta0block0', 'Delta1block0', 'Delta2block0']
    values = ['1600.0', '1.4', '0.1', '0.71', '0.003385', '500', '64', '64', '64', '2*M_PI/block0np0', '2*M_PI/block0np1', '2*M_PI/block0np2']
    simulation_name = 'opensbli'
    substitute_parameters(simulation_name, constants, values)

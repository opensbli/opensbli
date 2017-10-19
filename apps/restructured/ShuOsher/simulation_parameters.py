
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
    constants = ['gama', 'Minf', 'dt', 'niter', 'block0np0', 'Delta0block0']
    values = ['1.4', '0.1', '0.0002', 'ceil(1.8/dt)', '1600', '10.0/(block0np0-1)']
    simulation_name = 'opensbli'
    substitute_parameters(simulation_name, constants, values)

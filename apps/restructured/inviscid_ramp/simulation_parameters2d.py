
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
    constants = ['gama', 'Minf', 'dt', 'niter', 'block0np0', 'block0np1', 'Delta0block0', 'Delta1block0', 'lx', 'ly']
    values = ['1.4', '2.0', '0.005', '1000', '100', '100', '3.0/(block0np0-1)', '1.0/(block0np1-1)', "3.0", "1.0"]
    simulation_name = 'opensbli'
    substitute_parameters(simulation_name, constants, values)

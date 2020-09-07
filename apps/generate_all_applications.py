""" Script to generate all of the current OpenSBLI test cases."""
import os, subprocess
# List of the current applications
directories = [\
'/wave/',
'/euler_wave/',
'/shu_osher/',
'/Sod_shock_tube/',
'/taylor_green_vortex/',
'/taylor_green_vortex/TGsym/',
'/taylor_green_vortex/TGsym/',
'/viscous_shock_tube/',
'/kelvin_helmholtz/',
'/inviscid_shock_reflection/',
'/katzer_SBLI/',
'/channel_flow/laminar_2D/',
'/channel_flow/turbulent_3D/',
'/channel_flow/compressible_TCF_Central/',
'/channel_flow/compressible_TCF_TENO/',
'/Delery_bump/inviscid/',
'/Delery_bump/viscous/',
'/transitional_SBLI/',
'/cylinder/cylinder_central/',
'/cylinder/supersonic_cylinder/'
]
file_names = [\
'wave.py',
'euler_wave.py',
'shu_osher.py',
'Sod_shock_tube.py',
'taylor_green_vortex.py',
'TG_IsoT.py',
'TGsym.py',
'viscous_shock_tube.py',
'kelvin_helmholtz.py',
'inviscid_shock.py',
'katzer_SBLI.py',
'laminar_channel.py',
'turbulent_channel.py',
'turbulent_channel.py',
'turbulent_channel.py',
'inviscid_shock_delery_aerofoil_forced.py',
'viscous_shock_delery_aerofoil.py',
'transitional_SBLI.py',
'cylinder_central.py',
'supersonic_cylinder.py'
]

assert len(directories) == len(file_names)
print("Found %d OpenSBLI applications." % len(file_names))
# Current working directory
owd = os.getcwd()
# Optional diff between the generated codes
check_diff = False
if check_diff:
    import difflib
    # Set a directory containing previously generated C codes
    old_code_dir = os.environ['two'] + 'apps/'


with open(os.devnull, 'w') as devnull:

    for fname, directory in zip(file_names, directories):
        print("Generating the %s application." % (directory+fname))
        output_code = subprocess.call(["python %s" % fname], shell=True, cwd=owd+directory, stdout=devnull)
        if output_code == 0:
            print("%s generated successfully." % fname)
            # Compare the output code to a previously generated one
            if check_diff:
                file1, file2 = old_code_dir + directory + 'opensbli.cpp', owd + directory + 'opensbli.cpp'
                text1, text2 = open(file1).readlines(), open(file2).readlines()
                for line in difflib.unified_diff(text1, text2):
                    print(line)
        else:
            print("Generation of %s has failed." % fname)
            exit()

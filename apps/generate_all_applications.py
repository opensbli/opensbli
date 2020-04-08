""" Script to generate all of the current OpenSBLI test cases."""
import os, subprocess
# List of the current applications
directories = [\
'/wave/',
'/euler_wave/',
'/shu_osher/',
'/Sod_shock_tube/',
'/taylor_green_vortex/',
'/viscous_shock_tube/',
'/kelvin_helmholtz/',
'/inviscid_shock_reflection/',
'/katzer_SBLI/',
'/channel_flow/laminar_2D/',
'/channel_flow/turbulent_3D/',
'/channel_flow/compressible_TCF_Central/',
'/channel_flow/compressible_TCF_TENO/'
]
file_names = [\
'wave.py',
'euler_wave.py',
'shu_osher.py',
'Sod_shock_tube.py',
'taylor_green_vortex.py',
'viscous_shock_tube.py',
'kelvin_helmholtz.py',
'inviscid_shock.py',
'katzer_SBLI.py',
'laminar_channel.py',
'turbulent_channel.py',
'turbulent_channel.py',
'turbulent_channel.py'
]

assert len(directories) == len(file_names)
print("Found %d OpenSBLI applications." % len(file_names))
# Current working directory
owd = os.getcwd()

with open(os.devnull, 'w') as devnull:

    for fname, directory in zip(file_names, directories):
        print("Generating the %s application." % (directory+fname))
        output_code = subprocess.call(["python %s" % fname], shell=True, cwd=owd+directory, stdout=devnull)
        if output_code == 0:
            print("%s generated successfully." % fname)
        else:
            print("Generation of %s has failed." % fname)
            exit()

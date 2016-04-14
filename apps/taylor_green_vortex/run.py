import os, os.path
import sys
import subprocess
import shutil

base_name = "taylor_green_vortex"

# Generate the code
exit_code = subprocess.call("python %s.py" % base_name, shell=True)
if(exit_code != 0):
    print "Something went wrong when generating the code."
    sys.exit(1)

# Run the simulation
exit_code = subprocess.call("cd %s_opsc_code; make %s_cuda; ./%s_cuda" % (base_name, base_name, base_name), shell=True)

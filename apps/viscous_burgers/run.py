import os, os.path
import sys
import subprocess
import shutil

base_name = "viscous_burgers"

# Generate the code
exit_code = subprocess.call("python %s.py" % base_name, shell=True)
if(exit_code != 0):
    print "Something went wrong when generating the code."
    sys.exit(1)

# Run the simulation
exit_code = subprocess.call("cd %s_opsc_code; make %s_seq; ./%s_seq" % (base_name, base_name, base_name), shell=True)

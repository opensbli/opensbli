#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

import os, os.path
import sys
import subprocess
import shutil

base_name = "gaussian_bump"

# Generate the code
exit_code = subprocess.call("python %s.py" % base_name, shell=True)
if(exit_code != 0):
    print "Something went wrong when generating the code."
    sys.exit(1)

# Run the simulation
exit_code = subprocess.call("cd %s_opsc_code; make %s_seq; ./%s_seq" % (base_name, base_name, base_name), shell=True)

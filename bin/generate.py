#!/usr/bin/env python
import sys
import argparse

# Import local utility functions
from autofd.expansion import *

if(__name__ == "__main__"):
   # Parse the command line arguments provided by the user
   parser = argparse.ArgumentParser(prog="codegen", description="An automatic code generator which expands the equations written in Einstein notation, and writes out model code in OPSC or Fortran (serial) format.")
   parser.add_argument("equations_file", help="The path to the file containing the equations.", action="store", type=str)
   args = parser.parse_args()
   
   # Make this an absolute path
   abspath_to_equations_file = os.path.realpath(args.equations_file)
   
   # Expand out the equations
   expand_equations(abspath_to_equations_file)


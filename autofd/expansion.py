#!/usr/bin/env python
import sys
import datetime
import os
from re import *
import time

# Symbolic-related functions
from sympy import *
from sympy.parsing.sympy_parser import *
init_printing(use_latex=True)

# AutoFD functions
from .einstein_expansion import *
from .algorithm import *
from .utils import *

import logging
LOG = logging.getLogger(__name__)

BUILD_DIR = os.getcwd()

def read_input(equations_file):
    """ Read the equations file and extract the equations and their parameters. 
    
    :arg str equations_file: The path to the equations file.
    """

    # Remove leading and trailing white spaces and empty lines
    with open(equations_file) as f:
        read_file = [line for line in f.read().splitlines() if line]
    comment_lineno = [] # Line numbers of each (Python) comment line, which we want to ignore

    # Get all the comments in the file
    for ind,line in enumerate(read_file):
        if line[0] == '#':
            comment_lineno.append(ind)

    inp = EquationInput() # Define the inputs

    # Read inputs from the file
    inp.equations = read_file[comment_lineno[0]+1:comment_lineno[1]]
    inp.substitutions = read_file[comment_lineno[1]+1:comment_lineno[2]]
    inp.ndim = int(read_file[comment_lineno[2]+1])
    inp.constants = read_file[comment_lineno[3]+1:comment_lineno[4]]
    inp.coordinate_symbol = read_file[comment_lineno[4]+1:comment_lineno[5]]
    inp.metrics = read_file[comment_lineno[5]+1:comment_lineno[6]]
    inp.formulas = read_file[comment_lineno[6]+1:comment_lineno[7]]

    return inp

def expand_equations(equations_file):
   """ Perform an expansion of the equations, provided by the user, and written in Einstein notation.
   
   :arg str equations_file: The path to the equations file.
   """
   
   # Find out the path of the directory that the equations_file is in.
   base_path = os.path.dirname(equations_file)
   
   # Get the equation inputs
   inp = read_input(equations_file)

   # Find the tensor indices in the equations
   # Expand the equations
   start = time.time()
   equations = []
   for e in inp.equations:
     equations.append(Equation(e, inp))
   formulas = []
   for f in inp.formulas:
     formulas.append(Equation(f, inp))
   #for eq in eqs:
     #pprint(eq.expandedeq)
   end = time.time()
   LOG.debug('The time taken for tensor expansion of equations in %d Dimensions is %.2f seconds.' % (inp.ndim, end - start))
   
   # Prepare equations for algorithm
   algorithm_file_path = base_path + "/algorithm"
   read_file = [line for line in open(algorithm_file_path, "r").read().splitlines() if line]
   algorithm = read_alg(read_file)
   start = time.time()
   final_equation = PreparedEquations(equations, formulas, algorithm)
   end = time.time()
   LOG.debug('The time taken for tensor expansion of equations in %d Dimensions is %.2f seconds.' % (inp.ndim, end - start))

   # Output equations in LaTeX format.
   latex_file_path = BUILD_DIR + "/equations.tex"
   with open(latex_file_path, 'w') as latex_file:
      temp = []
      for e in equations:
        temp = temp + [e.parsed]
      inp = ['Equations', 'Satya P Jammy', 'University of Southampton']
      header, end = latex_article_header(inp)
      latex_file.write(header)
      temp = []
      for e in equations:
        temp = temp + [e.expandedeq]
      temp = flatten(temp)
      write_latex(latex_file, temp)
      latex_file.write(end)
   
   return

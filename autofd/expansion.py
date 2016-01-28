#!/usr/bin/env python
import os
from re import *
import time

# Symbolic-related functions
from sympy import *
from sympy.parsing.sympy_parser import *

# AutoFD functions
from .system import *
from .equations import *
from .algorithm import *
from .latex import LatexWriter

import logging
LOG = logging.getLogger(__name__)

BUILD_DIR = os.getcwd()


def expand_equations(equations_file):
    """ Perform an expansion of the equations, provided by the user, and written in Einstein notation.

    :arg str equations_file: The path to the equations file.
    """

    # Find out the path of the directory that the equations_file is in.
    base_path = os.path.dirname(equations_file)

    # Prepare the system and expand the equations
    system = System()
    system.read_input(equations_file)
    expanded_equations, expanded_formulas = system.expand_equations()

    # Output equations in LaTeX format.
    latex = LatexWriter()
    latex.open(path=BUILD_DIR + "/equations.tex")
    metadata = {"title": "Equations", "author": "Satya P Jammy", "institution": "University of Southampton"}
    latex.write_header(metadata)
    temp = flatten([e.expandedeq for e in expanded_equations])
    latex.write_equations(temp)
    latex.write_footer()
    latex.close()

    # Prepare the algorithm
    algorithm_file_path = base_path + "/algorithm"
    algorithm = Algorithm()
    algorithm.read_input(algorithm_file_path)

    # Generate the code
    start = time.time()
    final_equation = PreparedEquations(expanded_equations, expanded_formulas, algorithm)
    end = time.time()
    LOG.debug('The time taken to prepare the equations in %d Dimensions is %.2f seconds.' % (system.ndim, end - start))

    return

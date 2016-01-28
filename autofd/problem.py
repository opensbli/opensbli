#!/usr/bin/env python

import time

from .equations import *
from .algorithm import *

import logging
LOG = logging.getLogger(__name__)


class Problem(object):

    """ Describes the system we want to generate code for. """

    def __init__(self):
        self.equations = []
        self.substitutions = []
        self.ndim = []
        self.constants = []
        self.coordinate_symbol = []
        self.metrics = []
        self.formulas = []
        return

    def read_input(self, equation_file_path):
        """ Read the equations and algorithm files and extract their parameters.

        :arg str equation_file_path: The path to the equations file.
        """

        # Remove leading and trailing white spaces and empty lines
        with open(equation_file_path) as f:
            read_file = [line for line in f.read().splitlines() if line]
        comment_lineno = []  # Line numbers of each (Python) comment line, which we want to ignore

        # Get all the comments in the file
        for ind, line in enumerate(read_file):
            if line[0] == '#':
                comment_lineno.append(ind)

        # Read inputs from the file
        self.equations = read_file[comment_lineno[0]+1:comment_lineno[1]]
        self.substitutions = read_file[comment_lineno[1]+1:comment_lineno[2]]
        self.ndim = int(read_file[comment_lineno[2]+1])
        self.constants = read_file[comment_lineno[3]+1:comment_lineno[4]]
        self.coordinate_symbol = read_file[comment_lineno[4]+1:comment_lineno[5]]
        self.metrics = read_file[comment_lineno[5]+1:comment_lineno[6]]
        self.formulas = read_file[comment_lineno[6]+1:comment_lineno[7]]

        return

    def expand(self):
        """ Find the tensor indices in the equations, and then expand the equations. """

        LOG.info("Expanding equations...")

        start = time.time()
        expanded_equations = []
        for e in self.equations:
            expanded_equations.append(Equation(e, self))
        expanded_formulas = []
        for f in self.formulas:
            expanded_formulas.append(Equation(f, self))
        end = time.time()

        LOG.debug('The time taken for tensor expansion of equations in %d Dimensions is %.2f seconds.' % (self.ndim, end - start))
        return expanded_equations, expanded_formulas

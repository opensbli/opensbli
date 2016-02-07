#!/usr/bin/env python

import time

from .equations import *
from .algorithm import *

import logging
LOG = logging.getLogger(__name__)


class Problem(object):

    """ Describes the system we want to generate code for. """

    def __init__(self, equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas):
        self.equations = equations
        self.substitutions = substitutions
        self.ndim = ndim
        self.constants = constants
        self.coordinate_symbol = coordinate_symbol
        self.metrics = metrics
        self.formulas = formulas
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

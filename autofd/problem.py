#!/usr/bin/env python

#    AutoFD: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of AutoFD.

#    AutoFD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    AutoFD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with AutoFD.  If not, see <http://www.gnu.org/licenses/>.

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
            expanded_equations.append(Equation(e, self.ndim,self.coordinate_symbol, self.substitutions, self.constants))
        expanded_formulas = []
        for f in self.formulas:
            expanded_formulas.append(Equation(f, self.ndim, self.coordinate_symbol, self.substitutions, self.constants))
        end = time.time()

        LOG.debug('The time taken for tensor expansion of equations in %d Dimensions is %.2f seconds.' % (self.ndim, end - start))
        return expanded_equations, expanded_formulas

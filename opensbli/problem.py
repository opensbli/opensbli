#!/usr/bin/env python

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

import time

from .equations import *

import logging
LOG = logging.getLogger(__name__)


class Problem(object):

    """ Describes the system we want to generate code for. """

    def __init__(self, equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas):
        """ Store problem parameters, and create Equation objects for each user-provided equation/formula written in Einstien notation. """

        self.substitutions = substitutions
        self.ndim = ndim
        self.constants = constants
        self.coordinate_symbol = coordinate_symbol
        self.metrics = metrics

        LOG.info("Expanding equations...")
        start = time.time()
        # expand the equations
        self.equations = self.expand(equations)
        # expand the formulas
        self.formulas = self.expand(formulas)
        end = time.time()
        LOG.debug('The time taken for tensor expansion of equations in %d Dimensions is %.2f seconds.' % (self.ndim, end - start))
        return

    def expand(self, equations):
        """ Find the Einstein indices in the equations and formulas, and then expand them. """

        expanded = []
        for e in equations:
            expanded.append(Equation(e, self.ndim, self.coordinate_symbol, self.substitutions, self.constants))

        return expanded

    def get_expanded(self, equations):
        """ Return the lists of expanded equations and formulas. """

        expanded_equations = []
        for e in equations:
            expanded_equations.append(e.expanded)

        return expanded_equations

#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

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
        
        self.expand(equations, formulas)

        return

    def expand(self, equations, formulas):
        """ Find the Einstein indices in the equations and formulas, and then expand them. """

        LOG.info("Expanding equations...")

        start = time.time()

        self.equations = []
        for e in equations:
            self.equations.append(Equation(e, self.ndim, self.coordinate_symbol, self.substitutions, self.constants))
            
        self.formulas = []
        for f in formulas:
            self.formulas.append(Equation(f, self.ndim, self.coordinate_symbol, self.substitutions, self.constants))
            
        end = time.time()

        LOG.debug('The time taken for tensor expansion of equations in %d Dimensions is %.2f seconds.' % (self.ndim, end - start))
        return
        
    def get_expanded(self):  
        """ Return the lists of expanded equations and formulas. """
        
        expanded_equations = []
        for e in self.equations:
            expanded_equations.append(e.expanded)
        expanded_formulas = []
        for f in self.formulas:
            expanded_formulas.append(f.expanded)
        
        return expanded_equations, expanded_formulas

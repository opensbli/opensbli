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

from sympy import *
from sympy.parsing.sympy_parser import *

from .kernel import *
from .equations import *


class GridBasedInitialisation(object):

    """ Initialise the equations on the grid of solution points.  """

    def __init__(self, grid, ics):
        """ Create the initialisation kernels.

        :arg grid: The numerical grid of solution points.
        :arg list ics: A list of initial condition formulas.
        :returns: None
        """

        self.computations = []
        initialisation_equation = []
        for ic in ics:
            initialisation_equation.append(parse_expr(ic, local_dict={'grid': grid, 'Symbol': EinsteinTerm}))
        range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]

        self.computations.append(Kernel(initialisation_equation, range_of_evaluation, "Initialisation", grid))

        return

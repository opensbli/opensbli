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

import os
import subprocess
from sympy import *


class Evaluations(object):

    """ The evaluation of a LHS and RHS, containing information about what the LHS and RHS requires, whether there are any
    subevaluations that need doing (e.g. for second derivatives, e.g. d(du/dx)dy, du/dx should be evaluated first),
    and what the work arrays are. """

    def __init__(self, lhs, rhs, requires, subevals=None, wk=None):
        """ Set up the evaluation information. """

        if isinstance(lhs, Derivative):
            self.is_derivative = True
            self.is_formula = False

            if subevals:
                self.subevals = subevals
            else:
                self.subevals = [None]
            if wk:
                self.work = wk
            else:
                self.work = None

            self.formula = rhs
            self.requires = requires
            self.evaluation_range = []
        else:
            self.is_formula = True
            self.is_derivative = False

            if subevals:
                self.subevals = subevals
            else:
                self.subevals = [None]
            if wk:
                self.work = wk
            else:
                self.work = None

            self.formula = rhs
            self.requires = requires
            self.evaluation_range = []

        return

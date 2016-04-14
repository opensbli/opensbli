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


class FileIO(object):

    """ Saves the arrays provided after every n iterations. These will eventually be dumped into an HDF5 file. """

    def __init__(self, arrays, niter=None):
        """ Setup the 'save' arrays to dump to an HDF5 file.

        :arg arrays: The arrays to save to a file.
        :arg int niter: The number of iterations that should pass before the arrays are saved to a file. If niter is None, the arrays are saved at the end of the simulation.
        :returns: None
        """

        self.save_after = []
        self.save_arrays = []
        if niter == None:
            self.save_after += [True]
        else:
            # save after every n iterations and at the end of the simulation
            self.save_after += [niter, True]
        if isinstance(arrays, list):
            self.save_arrays += [arr.base for arr in arrays]
        else:
            self.save_arrays.append(arrays)
        return

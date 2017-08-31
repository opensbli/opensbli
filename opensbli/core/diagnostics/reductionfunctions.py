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

from opensbli.core.opensblifunctions import Function, BasicDiscretisation


class BlockReduction(BasicDiscretisation):
    pass


class BlockMax(Function, BasicDiscretisation):
    def __new__(cls):

        return


class BlockMin(Function, BasicDiscretisation):
    def __new__(cls):

        return


class BlockSum(Function, BasicDiscretisation):
    # This contains the sum of the variable
    pass


class BlockIntegral(Function, BasicDiscretisation):
    def __new__(cls, expr, *args):
        args = tuple(flatten([expr] + list(args)))
        # ret = super(CentralDerivative, cls).__new__(cls, *args, evaluate=False)
        ret.store = True  # By default all the derivatives are stored
        ret.local_evaluation = True
        return ret

        return

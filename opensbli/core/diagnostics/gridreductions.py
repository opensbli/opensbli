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

from opensbli.core.opensbliequations import NonSimulationEquations, Discretisation, Solution


class BlockReduction():
    def add_equations(cls, equation):
        equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            local = []
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq.lhs, eq.rhs)
                eq.set_vector(no)
                local += [eq]
            cls.equations += [local]
        else:
            equation = OpenSBLIEquation(equation.lhs, equation.rhs)
            cls.equations += [equation]
        return


class BlockMax(BlockReduction, Discretisation, Solution):
    def __new__(cls):

        return


class BlockMin(BlockReduction, Discretisation, Solution):
    def __new__(cls):

        return


class BlockSum(Discretisation, Solution):
    # This contains the sum of the variable
    pass


class BlockIntegral(BlockReduction, Discretisation, Solution):

    def __new__(cls):

        return

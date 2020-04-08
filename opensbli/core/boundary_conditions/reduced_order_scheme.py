from opensbli.core.opensbliobjects import DataSet
from sympy import Equality, Rational
from sympy.functions.elementary.piecewise import ExprCondPair


class ReducedOrder(object):
    """ Reduced order boundary scheme. Boundary point is 2nd order, one point above
    the wall is 3rd order."""

    def __init__(self):
        return

    def first_derivative(self, expression, direction, side):
        loc = list(list(expression.atoms(DataSet))[0].indices)
        weights = [[-Rational(3, 2), 2, -Rational(1, 2)], [-Rational(1, 3), -Rational(1, 2), 1, -Rational(1, 6)]]
        points = [[0, 1, 2], [-1, 0, 1, 2]]
        # weights = [[Rational(-11,6), 3, Rational(-3,2), Rational(1,3)], [-Rational(1,3), -Rational(1,2), 1, -Rational(1,6)]]
        # points = [[0, 1, 2, 3], [-1, 0, 1, 2]]
        if side == 0:
            fact = 1
        elif side == 1:
            fact = -1

        f_store = []
        for value in range(2):
            output_function = 0
            for i, index in enumerate(points[value]):
                new_loc = loc[:]
                new_loc[direction] += fact*index
                for dset in expression.atoms(DataSet):
                    expression = expression.replace(dset, dset.base[new_loc])
                output_function += expression*fact*weights[value][i]
            f_store.append(output_function)
        return f_store

    def second_derivative(self, expression, direction, side):
        loc = list(list(expression.atoms(DataSet))[0].indices)
        weights = [[2, -5, 4, -1], [Rational(11, 12), -Rational(5, 3), Rational(1, 2), Rational(1, 3), -Rational(1, 12)]]
        points = [[0, 1, 2, 3], [-1, 0, 1, 2, 3]]

        if side == 0:
            fact = 1
        elif side == 1:
            fact = -1

        f_store = []
        for value in range(2):
            output_function = 0
            for i, index in enumerate(points[value]):
                new_loc = loc[:]
                new_loc[direction] += fact*index
                for dset in expression.atoms(DataSet):
                    expression = expression.replace(dset, dset.base[new_loc])
                output_function += expression*fact*weights[value][i]
            f_store.append(output_function)
        return f_store

    def expr_cond_pairs(self, fn, direction, side, order, block):
        if order == 1:
            derivatives = self.first_derivative(fn, direction, side)
        elif order == 2:
            derivatives = self.second_derivative(fn, direction, side)
        else:
            raise NotImplementedError("Only first and second derivatives are implemented.")
        idx = block.grid_indexes[direction]
        if side == 0:
            mul_factor = 1
            start = block.ranges[direction][side]
        else:
            mul_factor = -1
            start = block.ranges[direction][side] - 1
        ecs = []
        for no, d in enumerate(derivatives):
            loc = start + mul_factor*no
            ecs += [ExprCondPair(d, Equality(idx, loc))]
        return ecs

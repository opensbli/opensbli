from opensbli.core.opensbliobjects import DataSet
from sympy import zeros, S, nsimplify, Matrix, Equality, Rational
from sympy.functions.elementary.piecewise import ExprCondPair


class ReducedAccess(object):
    """ Fourth order boundary scheme that uses a (-4,4) stencil."""

    def __init__(self):
        self.bc4_coefficients = self.first_der_coefficients()
        self.bc4_2_coefficients = self.second_der_coefficients()
        return

    def function_points(self, expression, direction, side):
        f_matrix = zeros(5, 5)
        loc = list(list(expression.atoms(DataSet))[0].indices)
        for shift in range(5):
            func_points = []
            for index in range(5):
                new_loc = loc[:]
                new_loc[direction] += index - shift
                for dset in expression.atoms(DataSet):
                    expression = expression.replace(dset, dset.base[new_loc])
                func_points.append(expression)
            f_matrix[:, shift] = func_points
        if side == 0:
            f_matrix = f_matrix[:, 0:2]
        elif side == 1:
            f_matrix = f_matrix.transpose()[:, 0:2]
        else:
            raise NotImplementedError("Side must be 0 or 1")
        return f_matrix

    def first_der_coefficients(self):
        coefficient_matrix = Matrix([[-25.0, 48.0, -36.0, 16.0, -3.0], [-3.0, -10.0, 18.0, -6.0, 1.0], [25.0, -48.0, 36.0, -16.0, 3.0],
                                     [3.0, 10.0, -18.0, 6.0, -1.0]])
        for i in range(coefficient_matrix.shape[0]):
            for j in range(coefficient_matrix.shape[1]):
                coefficient_matrix[i, j] = Rational(coefficient_matrix[i, j], 12.0)
        return coefficient_matrix

    def weight_function_points(self, func_points, direction, order, block, side, char_BC=False):
        if order == 1:
            if side == 0:
                coeffs = self.bc4_coefficients[0:2, :]
            elif side == 1:
                coeffs = self.bc4_coefficients[2:, :]
            weighted = zeros(2, 1)
            for i in range(2):
                weighted[i] = coeffs[i, :]*func_points[:, i]
        elif order == 2:
            h_sq = S.One  # The division of delta**2 is now applied in central derivative to reduce the divisions
            weighted = zeros(2, 1)
            for i in range(2):
                weighted[i] = h_sq*(self.bc4_2_coefficients[i, :]*func_points[0:5, i])
        else:
            raise NotImplementedError("Only 1st and 2nd derivatives implemented")
        return weighted

    def expr_cond_pairs(self, fn, direction, side, order, block):
        fn_pts = self.function_points(fn, direction, side)
        derivatives = self.weight_function_points(fn_pts, direction, order, block, side)
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

    def second_der_coefficients(self):
        """ Computes the finite-difference coefficients for the 2nd order one sided Carpenter wall boundary derivative.
        returns: Matrix: bc4_2: Matrix of stencil coefficients."""
        bc4_2 = Matrix([[35.0, -104.0, 114.0, -56.0, 11.0], [11.0, -20.0, 6.0, 4.0, -1.0]])/12.0
        for i in range(bc4_2.shape[0]):
            for j in range(bc4_2.shape[1]):
                bc4_2[i, j] = nsimplify(bc4_2[i, j])
        return bc4_2

from opensbli.core.opensbliobjects import DataSet
from sympy import zeros, S, nsimplify, sqrt, Matrix, Equality
from sympy.functions.elementary.piecewise import ExprCondPair


class Carpenter(object):
    """ 4th order one-sided Carpenter boundary treatment (https://doi.org/10.1006/jcph.1998.6114).
    If a boundary condition is an instance of ModifyCentralDerivative,
    central derivatives are replaced at that domain boundary by the Carpenter scheme."""

    def __init__(self):
        self.bc4_coefficients = self.carp4_coefficients()
        self.bc4_2_coefficients = self.second_der_coefficients()
        return

    def function_points(self, expression, direction, side):
        """ Create the function locations for evaluation of the Carpenter
        one-sided derivatives."""
        f_matrix = zeros(6, 6)
        loc = list(list(expression.atoms(DataSet))[0].indices)
        for shift in range(6):
            func_points = []
            for index in range(6):
                new_loc = loc[:]
                new_loc[direction] += index - shift
                for dset in expression.atoms(DataSet):
                    expression = expression.replace(dset, dset.base[new_loc])
                func_points.append(expression)
            f_matrix[:, shift] = func_points
        if side == 0:
            f_matrix = f_matrix[:, 0:4]
        elif side == 1:
            f_matrix = f_matrix.transpose()[:, 0:4]
        else:
            raise NotImplementedError("Side must be 0 or 1")
        return f_matrix

    def weight_function_points(self, func_points, direction, order, block, side):
        """ Multiply the function points with their weightings to form the derviatives."""
        if order == 1:
            h = S.One  # The division of delta is now applied in central derivative to reduce the divisions
            if side == 1:
                h = -S.One*h  # Modify the first derivatives for side ==1
            weighted = zeros(4, 1)
            for i in range(4):
                weighted[i] = h*(self.bc4_coefficients[i, :]*func_points[:, i])
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
        """ Computes the finite-difference coefficients for the 2nd derivative one sided Carpenter wall boundary closure.
        returns: Matrix: bc4_2: Matrix of stencil coefficients."""
        bc4_2 = Matrix([[35.0, -104.0, 114.0, -56.0, 11.0], [11.0, -20.0, 6.0, 4.0, -1.0]])/12.0
        for i in range(bc4_2.shape[0]):
            for j in range(bc4_2.shape[1]):
                bc4_2[i, j] = nsimplify(bc4_2[i, j])
        return bc4_2

    def carp4_coefficients(self):
        """ Computes the finite-difference coefficients for the 1st order one sided Carpenter wall boundary derivative.
        :returns: Matrix: bc4: Matrix of stencil coefficients."""
        R1 = -(2177.0*sqrt(295369.0)-1166427.0)/25488.0
        R2 = (66195.0*sqrt(53.0)*sqrt(5573.0)-35909375.0)/101952.0

        al4_0 = [-(216.0*R2+2160.0*R1-2125.0)/12960.0, (81.0*R2+675.0*R1+415.0)/540.0, -(72.0*R2+720.0*R1+445.0)/1440.0, -(108.0*R2+756.0*R1+421.0)/1296.0]
        al4_1 = [(81.0*R2+675.0*R1+415.0)/540.0, -(4104.0*R2+32400.0*R1+11225.0)/4320.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(216.0*R2+2160.0*R1+655.0)/4320.0]
        al4_2 = [-(72.0*R2+720.0*R1+445.0)/1440.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(4104.0*R2+32400.0*R1+12785.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0]
        al4_3 = [-(108.0*R2+756.0*R1+421.0)/1296.0, -(216.0*R2+2160.0*R1+655.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0, -(216.0*R2+2160.0*R1-12085.0)/12960.0]

        al4 = Matrix([al4_0, al4_1, al4_2, al4_3])

        ar4_0 = [(-1.0)/2.0, -(864.0*R2+6480.0*R1+305.0)/4320.0, (216.0*R2+1620.0*R1+725.0)/540.0, -(864.0*R2+6480.0*R1+3335.0)/4320.0, 0.0, 0.0]
        ar4_1 = [(864.0*R2+6480.0*R1+305.0)/4320.0, 0.0, -(864.0*R2+6480.0*R1+2315.0)/1440.0, (108.0*R2+810.0*R1+415.0)/270.0, 0.0, 0.0]
        ar4_2 = [-(216.0*R2+1620.0*R1+725.0)/540.0, (864.0*R2+6480.0*R1+2315.0)/1440.0, 0.0, -(864.0*R2+6480.0*R1+785.0)/4320.0, -1.0/12.0, 0.0]
        ar4_3 = [(864.0*R2+6480.0*R1+3335.0)/4320.0, -(108.0*R2+810.0*R1+415.0)/270.0, (864.0*R2+6480.0*R1+785.0)/4320.0, 0.0, 8.0/12.0, -1.0/12.0]
        ar4 = Matrix([ar4_0, ar4_1, ar4_2, ar4_3])
        # Form inverse and convert to rational
        al4_inv = al4.inv()
        bc4 = al4_inv*ar4
        return bc4

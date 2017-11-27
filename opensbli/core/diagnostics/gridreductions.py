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

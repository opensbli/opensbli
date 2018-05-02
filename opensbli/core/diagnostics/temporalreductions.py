from opensbli.core.opensbliequations import NonSimulationEquations, Discretisation, Solution


class BlockReduction(object):
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


class TemporalMax(BlockReduction, Discretisation, Solution):
    def __new__(cls):

        return


class TemporalMin(BlockReduction, Discretisation, Solution):
    def __new__(cls):

        return


class TemporalSum(Discretisation, Solution, BlockReduction):
    # This contains the sum of the variable
    pass


class TemporalMean(Discretisation, Solution):

    def __new__(cls):
        # TempralMean(a) --> Mean(OpenSBLISum(a, temporal), temoral)
        # Eq(meanrhou_i_j, TempralMean(rhou_i*rhou_j/rho))
        # [Eq(meanrhou0,  OpenSBLIMean(OpenSBLISum(CD(rhou0,x0), temporal), temporal))]
        # Eq(meanrhou_i, SpatialMean(rhou_i))
        # [Eq(meanrhou0,  OpenSBLIMean(OpenSBLISum(rhou0, spatial), spatial))]

        return
    
    def kernels()

a = Eq(mu0, TempralMean(u0))

diagn = NonSimulationEquations()

diagn.add_equatiosn(a)

Eq(Mean(mu0, temporal), Mean(Sum(u0, temporal), temproal))

input = "Eq(mu0, Mean(u0, temporal))"
#the above input would be transformed to this while parsing
Eq(DataObject("mu0"), TemporalMean(TemporalS(DataObject("u0"))))

input = "Eq(mu0, Mean(u0))" # This would be spatial mean
# It would be converted to 
Eq(ReductionVariable("mu0"), SpatianMean(SpatialSum(DataObject("u0"))))

Eq(sumu0,  OpenSBLISum(u0, spatial))
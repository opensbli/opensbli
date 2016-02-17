#!/usr/bin/env python

import os
import pytest

from autofd.equations import *


@pytest.fixture
def mass():
    return Equation("Eq(Der(rho,t), -Conservative(rhou_i,x_j))", 3, substitutions=[], constants=[])


@pytest.fixture
def stress_tensor():
    return "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j) + Der(u_j,x_i) - (2/3) * KroneckerDelta(_i,_j)*Der(u_k,x_k)))"


@pytest.fixture
def heat_flux():
    return "Eq(q_i, -(mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_i))"


@pytest.fixture
def momentum(stress_tensor):
    return Equation("Eq(Der(rhou_i,t) ,-Conservative(rhou_i*u_j + p* KroneckerDelta(_i,_j),x_j) + Der(tau_i_j,x_j) )", 3, substitutions=[stress_tensor], constants=["Re", "mu"])


@pytest.fixture
def energy(stress_tensor, heat_flux):
    return Equation("Eq(Der(rhoE,t),-Conservative((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )", 3, substitutions=[stress_tensor, heat_flux], constants=["Re", "Pr", "mu", "gama", "Minf"])


@pytest.fixture
def c():
    """ A scalar EinsteinTerm """
    return EinsteinTerm("c")


@pytest.fixture
def c_constant():
    """ A constant scalar EinsteinTerm """
    c = EinsteinTerm("c")
    c.constant = True
    return c


@pytest.fixture
def u_i():
    """ A vector EinsteinTerm """
    return EinsteinTerm("u_i")


@pytest.fixture
def tau_i_j():
    """ A tensor EinsteinTerm """
    return EinsteinTerm("tau_i_j")


def test_indices(c, u_i, tau_i_j):
    """ Check that all indices of an Einstein term are parsed correctly. """

    assert c.indices == []  # Scalar
    assert u_i.indices == ["_i"]  # Vector
    assert tau_i_j.indices == ["_i", "_j"]  # Tensor


def test_name(c, u_i, tau_i_j):
    """ Check that the name of an Einstein term is correct. """

    assert c.name == "c"  # Scalar
    assert u_i.name == "u_i"  # Vector
    assert tau_i_j.name == "tau_i_j"  # Tensor


def test_has_index(c, u_i, tau_i_j):
    """ Check that various Einstein terms have the indices that they should have. """

    assert not c.has_index("_l")

    assert u_i.has_index("_i")
    assert not u_i.has_index("_j")

    assert tau_i_j.has_index("_i")
    assert tau_i_j.has_index("_j")
    assert not tau_i_j.has_index("_k")


def test_expand_index(c_constant, u_i, tau_i_j):
    """ Check that an index can successfully be applied to an Einstein term. """

    new = c_constant.expand_index("_i", 0)
    assert new == c_constant  # This should stay the same since the Einstein term is constant.

    new = u_i.expand_index("_i", 0)
    assert not new.has_index("_i")

    new = tau_i_j.expand_index("_i", 0)
    assert new.has_index("_j")
    assert not new.has_index("_i")


def test_equation_expand(mass, momentum, energy):
    """ Check that the Navier-Stokes equations are expanded correctly. """

    # Conservation of mass
    assert str(mass.parsed) == "Der(rho, t) == -Conservative(rhou_i, x_j)"
    assert isinstance(mass.expanded, list)
    assert len(mass.expanded) == 1
    assert str(mass.expanded[0]) == "Derivative(rho(t, x0, x1, x2), t) == -Derivative(rhou0(t, x0, x1, x2), x0) - Derivative(rhou0(t, x0, x1, x2), x1) - Derivative(rhou0(t, x0, x1, x2), x2) - Derivative(rhou1(t, x0, x1, x2), x0) - Derivative(rhou1(t, x0, x1, x2), x1) - Derivative(rhou1(t, x0, x1, x2), x2) - Derivative(rhou2(t, x0, x1, x2), x0) - Derivative(rhou2(t, x0, x1, x2), x1) - Derivative(rhou2(t, x0, x1, x2), x2)"
    assert str(mass.expanded[0]).count("Derivative") == 10  # 1 time Derivative, then 3 Derivative for each of the 3 dimensions

    # Conservation of momentum
    assert str(momentum.parsed) == "Der(rhou_i, t) == -Conservative(p*KroneckerDelta(_i, _j) + rhou_i*u_j, x_j) + Der(mu*(Der(u_i, x_j) + Der(u_j, x_i) - 2*Der(u_k, x_k)*KroneckerDelta(_i, _j)/3)/Re, x_j)"
    assert isinstance(momentum.expanded, list)
    assert len(momentum.expanded) == 3  # Velocity is a vector, and the problem is 3D, so we should get 3 equations back (one for each component of the vector).
    assert str(momentum.expanded[0]) == "Derivative(rhou0(t, x0, x1, x2), t) == -Derivative(rhou0(t, x0, x1, x2)*u1(t, x0, x1, x2), x1) - Derivative(rhou0(t, x0, x1, x2)*u2(t, x0, x1, x2), x2) - Derivative(p(t, x0, x1, x2) + rhou0(t, x0, x1, x2)*u0(t, x0, x1, x2), x0) + mu*(Derivative(u0(t, x0, x1, x2), x1, x1) + Derivative(u1(t, x0, x1, x2), x0, x1))/Re + mu*(Derivative(u0(t, x0, x1, x2), x2, x2) + Derivative(u2(t, x0, x1, x2), x0, x2))/Re + mu*(4*Derivative(u0(t, x0, x1, x2), x0, x0)/3 - 2*Derivative(u1(t, x0, x1, x2), x0, x1)/3 - 2*Derivative(u2(t, x0, x1, x2), x0, x2)/3)/Re"
    assert str(momentum.expanded[1]) == "Derivative(rhou1(t, x0, x1, x2), t) == -Derivative(rhou1(t, x0, x1, x2)*u0(t, x0, x1, x2), x0) - Derivative(rhou1(t, x0, x1, x2)*u2(t, x0, x1, x2), x2) - Derivative(p(t, x0, x1, x2) + rhou1(t, x0, x1, x2)*u1(t, x0, x1, x2), x1) + mu*(Derivative(u0(t, x0, x1, x2), x0, x1) + Derivative(u1(t, x0, x1, x2), x0, x0))/Re + mu*(Derivative(u1(t, x0, x1, x2), x2, x2) + Derivative(u2(t, x0, x1, x2), x1, x2))/Re + mu*(-2*Derivative(u0(t, x0, x1, x2), x0, x1)/3 + 4*Derivative(u1(t, x0, x1, x2), x1, x1)/3 - 2*Derivative(u2(t, x0, x1, x2), x1, x2)/3)/Re"
    assert str(momentum.expanded[2]) == "Derivative(rhou2(t, x0, x1, x2), t) == -Derivative(rhou2(t, x0, x1, x2)*u0(t, x0, x1, x2), x0) - Derivative(rhou2(t, x0, x1, x2)*u1(t, x0, x1, x2), x1) - Derivative(p(t, x0, x1, x2) + rhou2(t, x0, x1, x2)*u2(t, x0, x1, x2), x2) + mu*(Derivative(u0(t, x0, x1, x2), x0, x2) + Derivative(u2(t, x0, x1, x2), x0, x0))/Re + mu*(Derivative(u1(t, x0, x1, x2), x1, x2) + Derivative(u2(t, x0, x1, x2), x1, x1))/Re + mu*(-2*Derivative(u0(t, x0, x1, x2), x0, x2)/3 - 2*Derivative(u1(t, x0, x1, x2), x1, x2)/3 + 4*Derivative(u2(t, x0, x1, x2), x2, x2)/3)/Re"

    # Conservation of energy
    print energy.parsed
    print str(energy.expanded[0]) 
    assert str(energy.parsed) == "Der(rhoE, t) == -Conservative(u_j*(p + rhoE), x_j) + Der(mu*u_i*(Der(u_i, x_j) + Der(u_j, x_i) - 2*Der(u_k, x_k)*KroneckerDelta(_i, _j)/3)/Re, x_j) - Der(-mu*Der(T, x_i)/(Minf**2*Pr*Re*(gama - 1)), x_i)"
    assert isinstance(energy.expanded, list)
    assert len(energy.expanded) == 1
    assert str(energy.expanded[0]) == "Derivative(rhoE(t, x0, x1, x2), t) == -Derivative((p(t, x0, x1, x2) + rhoE(t, x0, x1, x2))*u0(t, x0, x1, x2), x0) - Derivative((p(t, x0, x1, x2) + rhoE(t, x0, x1, x2))*u1(t, x0, x1, x2), x1) - Derivative((p(t, x0, x1, x2) + rhoE(t, x0, x1, x2))*u2(t, x0, x1, x2), x2) + mu*((Derivative(u0(t, x0, x1, x2), x1) + Derivative(u1(t, x0, x1, x2), x0))*Derivative(u0(t, x0, x1, x2), x1) + (Derivative(u1(t, x0, x1, x2), x2) + Derivative(u2(t, x0, x1, x2), x1))*Derivative(u2(t, x0, x1, x2), x1) + (Derivative(u0(t, x0, x1, x2), x1, x1) + Derivative(u1(t, x0, x1, x2), x0, x1))*u0(t, x0, x1, x2) + (Derivative(u1(t, x0, x1, x2), x1, x2) + Derivative(u2(t, x0, x1, x2), x1, x1))*u2(t, x0, x1, x2) + (-2*Derivative(u0(t, x0, x1, x2), x0)/3 + 4*Derivative(u1(t, x0, x1, x2), x1)/3 - 2*Derivative(u2(t, x0, x1, x2), x2)/3)*Derivative(u1(t, x0, x1, x2), x1) + (-2*Derivative(u0(t, x0, x1, x2), x0, x1)/3 + 4*Derivative(u1(t, x0, x1, x2), x1, x1)/3 - 2*Derivative(u2(t, x0, x1, x2), x1, x2)/3)*u1(t, x0, x1, x2))/Re + mu*((Derivative(u0(t, x0, x1, x2), x1) + Derivative(u1(t, x0, x1, x2), x0))*Derivative(u1(t, x0, x1, x2), x0) + (Derivative(u0(t, x0, x1, x2), x2) + Derivative(u2(t, x0, x1, x2), x0))*Derivative(u2(t, x0, x1, x2), x0) + (Derivative(u0(t, x0, x1, x2), x0, x1) + Derivative(u1(t, x0, x1, x2), x0, x0))*u1(t, x0, x1, x2) + (Derivative(u0(t, x0, x1, x2), x0, x2) + Derivative(u2(t, x0, x1, x2), x0, x0))*u2(t, x0, x1, x2) + (4*Derivative(u0(t, x0, x1, x2), x0)/3 - 2*Derivative(u1(t, x0, x1, x2), x1)/3 - 2*Derivative(u2(t, x0, x1, x2), x2)/3)*Derivative(u0(t, x0, x1, x2), x0) + (4*Derivative(u0(t, x0, x1, x2), x0, x0)/3 - 2*Derivative(u1(t, x0, x1, x2), x0, x1)/3 - 2*Derivative(u2(t, x0, x1, x2), x0, x2)/3)*u0(t, x0, x1, x2))/Re + mu*((Derivative(u0(t, x0, x1, x2), x2) + Derivative(u2(t, x0, x1, x2), x0))*Derivative(u0(t, x0, x1, x2), x2) + (Derivative(u1(t, x0, x1, x2), x2) + Derivative(u2(t, x0, x1, x2), x1))*Derivative(u1(t, x0, x1, x2), x2) + (Derivative(u0(t, x0, x1, x2), x2, x2) + Derivative(u2(t, x0, x1, x2), x0, x2))*u0(t, x0, x1, x2) + (Derivative(u1(t, x0, x1, x2), x2, x2) + Derivative(u2(t, x0, x1, x2), x1, x2))*u1(t, x0, x1, x2) + (-2*Derivative(u0(t, x0, x1, x2), x0)/3 - 2*Derivative(u1(t, x0, x1, x2), x1)/3 + 4*Derivative(u2(t, x0, x1, x2), x2)/3)*Derivative(u2(t, x0, x1, x2), x2) + (-2*Derivative(u0(t, x0, x1, x2), x0, x2)/3 - 2*Derivative(u1(t, x0, x1, x2), x1, x2)/3 + 4*Derivative(u2(t, x0, x1, x2), x2, x2)/3)*u2(t, x0, x1, x2))/Re + mu*Derivative(T(t, x0, x1, x2), x0, x0)/(Minf**2*Pr*Re*(gama - 1)) + mu*Derivative(T(t, x0, x1, x2), x1, x1)/(Minf**2*Pr*Re*(gama - 1)) + mu*Derivative(T(t, x0, x1, x2), x2, x2)/(Minf**2*Pr*Re*(gama - 1))"


def test_has_nested_derivative():
    """ Ensure that nested derivatives (Der) objects are correctly found. """

    ndim = 2
    equation = Equation("Eq(Der(rhoE,t), Der(u_i*Der(u_i,x_j),x_j) )", ndim)
    expansion = EinsteinExpansion(equation.parsed, ndim)
    assert expansion.has_nested_derivative(equation.parsed.rhs)
    assert not expansion.has_nested_derivative(equation.parsed.lhs)
    

def test_get_nested():
    """ Ensure that all terms (whole tree nodes) with index "i" are identified correctly. """
    ndim = 2
    equation = Equation("Eq(Der(rhoE,t), Der(u_i*Der(u_i,x_j),x_j) )", ndim)
    expansion = EinsteinExpansion(equation.parsed, ndim)
    nested_terms = expansion.get_nested(equation.parsed.rhs, "_i")
    assert len(nested_terms) == 1
    assert str(nested_terms[0]) == "u_i*Der(u_i, x_j)"


def test_get_einstein_indices():
    """ Get all the Einstein indices in the equation's RHS. """
    ndim = 2
    equation = Equation("Eq(Der(rhoE,t), Der(u_i*Der(u_i,x_j),x_j) )", ndim)
    expansion = EinsteinExpansion(equation.parsed, ndim)
    indices = expansion.get_einstein_indices(expression = equation.parsed.rhs)
    assert len(indices) == 2
    assert indices == set(["_i", "_j"])


def test_equations_to_dict(mass):
    """ Ensure that a list of SymPy Eq objects is converted to a dictionary. """

    equations = [mass.expanded[0]]
    assert isinstance(equations, list)
    equations = equations_to_dict(equations)
    assert isinstance(equations, dict)


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))

#!/usr/bin/env python

import os
import pytest

from sympy import Symbol, Idx

from autofd.equations import *

@pytest.fixture
def coordinate_symbol():
    return "x"

@pytest.fixture
def mass(coordinate_symbol):
    return Equation("Eq(Der(rho,t), -Conservative(rhou_j,x_j))", 3, coordinate_symbol, substitutions=[], constants=[])


@pytest.fixture
def stress_tensor():
    return "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j) + Der(u_j,x_i) - (2/3) * KD(_i,_j)*Der(u_k,x_k)))"


@pytest.fixture
def heat_flux():
    return "Eq(q_i, -(mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_i))"


@pytest.fixture
def momentum(stress_tensor, coordinate_symbol):
    return Equation("Eq(Der(rhou_i,t) ,-Conservative(rhou_i*u_j + p* KD(_i,_j),x_j) + Der(tau_i_j,x_j) )", 3, coordinate_symbol, substitutions=[stress_tensor], constants=["Re", "mu"])


@pytest.fixture
def energy(stress_tensor, heat_flux, coordinate_symbol):
    return Equation("Eq(Der(rhoE,t),-Conservative((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )", 3, coordinate_symbol, substitutions=[stress_tensor, heat_flux], constants=["Re", "Pr", "mu", "gama", "Minf"])


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


def test_get_indices(c, u_i, tau_i_j):
    """ Check that all indices of an EinsteinTerm are parsed correctly. """

    assert c.get_indices() == []  # Scalar
    assert u_i.get_indices() == [Idx(Symbol("i", integer=True))]  # Vector
    assert tau_i_j.get_indices() == [Idx(Symbol("i", integer=True)), Idx(Symbol("j", integer=True))]  # Tensor


def test_name(c, u_i, tau_i_j):
    """ Check that the base name of an EinsteinTerm is correct. """

    assert c.get_base() == "c"  # Scalar
    assert u_i.get_base() == "u"  # Vector
    assert tau_i_j.get_base() == "tau"  # Tensor


def test_equation_expand(mass, momentum, energy):
    """ Check that the Navier-Stokes equations are expanded correctly. """

    # Conservation of mass
    assert str(mass.parsed) == "Der(rho, t) == -Conservative(rhou_j, x_j)"
    assert isinstance(mass.expanded, list)
    assert len(mass.expanded) == 1
    assert str(mass.expanded[0]) == "Derivative(rho[x0, x1, x2, t], t) == -Derivative(rhou0[x0, x1, x2, t], x0) - Derivative(rhou1[x0, x1, x2, t], x1) - Derivative(rhou2[x0, x1, x2, t], x2)"
    assert str(mass.expanded[0]).count("Derivative") == 4  # 1 time Derivative, then 3 Derivative for each dimension



#def test_get_einstein_indices():
#    """ Get all the Einstein indices in the equation's RHS. """
#    ndim = 2
#    equation = Equation("Eq(Der(rhoE,t), Der(u_i*Der(u_i,x_j),x_j) )", ndim)
#    expansion = EinsteinExpansion(equation.parsed, ndim)
#    indices = expansion.get_einstein_indices(expression=equation.parsed.rhs)
#    assert len(indices) == 2
#    assert indices == set(["_i", "_j"])


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))

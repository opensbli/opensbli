#!/usr/bin/env python

import os
import pytest

from sympy import Symbol, Idx


# OpenSBLI functions
from opensbli.equations import *

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
    return Equation("Eq(Der(rhoE,t),-Conservative((p+rhoE)*u_j,x_j) -Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )", 3, coordinate_symbol, substitutions=[stress_tensor, heat_flux], constants=["Re", "Pr", "mu", "gama", "Minf"])


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


def test_get_base(c, u_i, tau_i_j):
    """ Check that the base name of an EinsteinTerm is correct. """

    assert c.get_base() == "c"  # Scalar
    assert u_i.get_base() == "u"  # Vector
    assert tau_i_j.get_base() == "tau"  # Tensor


def test_get_expanded(c, u_i, tau_i_j):
    """ Check that the indices of an EinsteinTerm are expanded correctly to produce a new EinsteinTerm object. """

    index_map = [(Idx(Symbol("i", integer=True)), 0), (Idx(Symbol("j", integer=True)), 1)]
    assert c.get_expanded(index_map) == c  # Scalar. No indices should be expanded here.
    assert u_i.get_expanded(index_map) == EinsteinTerm("u0")  # Vector.
    assert tau_i_j.get_expanded(index_map) == EinsteinTerm("tau01")  # Tensor.


def test_equation(mass, momentum, energy):
    """ Check that the Navier-Stokes equations are expanded correctly. """

    # Conservation of mass
    assert str(mass.parsed.lhs) == "Der(rho, t)"
    assert str(mass.parsed.rhs)== "-Conservative(rhou_j, x_j)"
    assert isinstance(mass.expanded, list)
    assert len(mass.expanded) == 1
    assert str(mass.expanded[0].lhs) == "Derivative(rho[x0, x1, x2, t], t)"
    assert str(mass.expanded[0].rhs) == "-Derivative(rhou0[x0, x1, x2, t], x0) - Derivative(rhou1[x0, x1, x2, t], x1) - Derivative(rhou2[x0, x1, x2, t], x2)"
    assert str(mass.expanded[0]).count("Derivative") == 4  # 1 time Derivative, then 3 Derivative for each dimension

    # Conservation of momentum
    assert str(momentum.parsed.lhs) == "Der(rhou_i, t)"
    assert str(momentum.parsed.rhs) == "-Conservative(p*KD(_i, _j) + rhou_i*u_j, x_j) + Der(mu*Re**(-1)*(Der(u_i, x_j) + Der(u_j, x_i) - 2*KD(_i, _j)*Der(u_k, x_k)/3), x_j)"
    assert isinstance(momentum.expanded, list)
    assert len(momentum.expanded) == 3  # Velocity is a vector, and the problem is 3D, so we should get 3 equations back (one for each component of the vector).
    assert str(momentum.expanded[0].lhs) == "Derivative(rhou0[x0, x1, x2, t], t)"
    assert str(momentum.expanded[0].rhs) == "mu*(Derivative(u0[x0, x1, x2, t], x1, x1) + Derivative(u1[x0, x1, x2, t], x0, x1))*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x2, x2) + Derivative(u2[x0, x1, x2, t], x0, x2))*Re**(-1) + mu*(4*Derivative(u0[x0, x1, x2, t], x0, x0)/3 - 2*Derivative(u1[x0, x1, x2, t], x0, x1)/3 - 2*Derivative(u2[x0, x1, x2, t], x0, x2)/3)*Re**(-1) - Derivative(rhou0[x0, x1, x2, t]*u1[x0, x1, x2, t], x1) - Derivative(rhou0[x0, x1, x2, t]*u2[x0, x1, x2, t], x2) - Derivative(p[x0, x1, x2, t] + rhou0[x0, x1, x2, t]*u0[x0, x1, x2, t], x0)"
    assert str(momentum.expanded[1].lhs) == "Derivative(rhou1[x0, x1, x2, t], t)"
    assert str(momentum.expanded[1].rhs) == "mu*(Derivative(u0[x0, x1, x2, t], x0, x1) + Derivative(u1[x0, x1, x2, t], x0, x0))*Re**(-1) + mu*(Derivative(u1[x0, x1, x2, t], x2, x2) + Derivative(u2[x0, x1, x2, t], x1, x2))*Re**(-1) + mu*(-2*Derivative(u0[x0, x1, x2, t], x0, x1)/3 + 4*Derivative(u1[x0, x1, x2, t], x1, x1)/3 - 2*Derivative(u2[x0, x1, x2, t], x1, x2)/3)*Re**(-1) - Derivative(rhou1[x0, x1, x2, t]*u0[x0, x1, x2, t], x0) - Derivative(rhou1[x0, x1, x2, t]*u2[x0, x1, x2, t], x2) - Derivative(p[x0, x1, x2, t] + rhou1[x0, x1, x2, t]*u1[x0, x1, x2, t], x1)"
    assert str(momentum.expanded[2].lhs) == "Derivative(rhou2[x0, x1, x2, t], t)"
    assert str(momentum.expanded[2].rhs) == "mu*(Derivative(u0[x0, x1, x2, t], x0, x2) + Derivative(u2[x0, x1, x2, t], x0, x0))*Re**(-1) + mu*(Derivative(u1[x0, x1, x2, t], x1, x2) + Derivative(u2[x0, x1, x2, t], x1, x1))*Re**(-1) + mu*(-2*Derivative(u0[x0, x1, x2, t], x0, x2)/3 - 2*Derivative(u1[x0, x1, x2, t], x1, x2)/3 + 4*Derivative(u2[x0, x1, x2, t], x2, x2)/3)*Re**(-1) - Derivative(rhou2[x0, x1, x2, t]*u0[x0, x1, x2, t], x0) - Derivative(rhou2[x0, x1, x2, t]*u1[x0, x1, x2, t], x1) - Derivative(p[x0, x1, x2, t] + rhou2[x0, x1, x2, t]*u2[x0, x1, x2, t], x2)"

    # Conservation of energy
    assert str(energy.parsed.lhs) == "Der(rhoE, t)"
    assert str(energy.parsed.rhs) == "-Conservative((p + rhoE)*u_j, x_j) - Der(q_j, x_j) + Der(mu*u_i*Re**(-1)*(Der(u_i, x_j) + Der(u_j, x_i) - 2*KD(_i, _j)*Der(u_k, x_k)/3), x_j)"
    assert isinstance(energy.expanded, list)
    assert len(energy.expanded) == 1
    assert str(energy.expanded[0].lhs) == "Derivative(rhoE[x0, x1, x2, t], t)"
    assert str(energy.expanded[0].rhs) == "mu*(Derivative(u0[x0, x1, x2, t], x1) + Derivative(u1[x0, x1, x2, t], x0))*Derivative(u0[x0, x1, x2, t], x1)*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x1) + Derivative(u1[x0, x1, x2, t], x0))*Derivative(u1[x0, x1, x2, t], x0)*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x2) + Derivative(u2[x0, x1, x2, t], x0))*Derivative(u0[x0, x1, x2, t], x2)*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x2) + Derivative(u2[x0, x1, x2, t], x0))*Derivative(u2[x0, x1, x2, t], x0)*Re**(-1) + mu*(Derivative(u1[x0, x1, x2, t], x2) + Derivative(u2[x0, x1, x2, t], x1))*Derivative(u1[x0, x1, x2, t], x2)*Re**(-1) + mu*(Derivative(u1[x0, x1, x2, t], x2) + Derivative(u2[x0, x1, x2, t], x1))*Derivative(u2[x0, x1, x2, t], x1)*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x0, x1) + Derivative(u1[x0, x1, x2, t], x0, x0))*u1[x0, x1, x2, t]*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x0, x2) + Derivative(u2[x0, x1, x2, t], x0, x0))*u2[x0, x1, x2, t]*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x1, x1) + Derivative(u1[x0, x1, x2, t], x0, x1))*u0[x0, x1, x2, t]*Re**(-1) + mu*(Derivative(u0[x0, x1, x2, t], x2, x2) + Derivative(u2[x0, x1, x2, t], x0, x2))*u0[x0, x1, x2, t]*Re**(-1) + mu*(Derivative(u1[x0, x1, x2, t], x1, x2) + Derivative(u2[x0, x1, x2, t], x1, x1))*u2[x0, x1, x2, t]*Re**(-1) + mu*(Derivative(u1[x0, x1, x2, t], x2, x2) + Derivative(u2[x0, x1, x2, t], x1, x2))*u1[x0, x1, x2, t]*Re**(-1) + mu*(-2*Derivative(u0[x0, x1, x2, t], x0)/3 - 2*Derivative(u1[x0, x1, x2, t], x1)/3 + 4*Derivative(u2[x0, x1, x2, t], x2)/3)*Derivative(u2[x0, x1, x2, t], x2)*Re**(-1) + mu*(-2*Derivative(u0[x0, x1, x2, t], x0)/3 + 4*Derivative(u1[x0, x1, x2, t], x1)/3 - 2*Derivative(u2[x0, x1, x2, t], x2)/3)*Derivative(u1[x0, x1, x2, t], x1)*Re**(-1) + mu*(4*Derivative(u0[x0, x1, x2, t], x0)/3 - 2*Derivative(u1[x0, x1, x2, t], x1)/3 - 2*Derivative(u2[x0, x1, x2, t], x2)/3)*Derivative(u0[x0, x1, x2, t], x0)*Re**(-1) + mu*(4*Derivative(u0[x0, x1, x2, t], x0, x0)/3 - 2*Derivative(u1[x0, x1, x2, t], x0, x1)/3 - 2*Derivative(u2[x0, x1, x2, t], x0, x2)/3)*u0[x0, x1, x2, t]*Re**(-1) + mu*(-2*Derivative(u0[x0, x1, x2, t], x0, x1)/3 + 4*Derivative(u1[x0, x1, x2, t], x1, x1)/3 - 2*Derivative(u2[x0, x1, x2, t], x1, x2)/3)*u1[x0, x1, x2, t]*Re**(-1) + mu*(-2*Derivative(u0[x0, x1, x2, t], x0, x2)/3 - 2*Derivative(u1[x0, x1, x2, t], x1, x2)/3 + 4*Derivative(u2[x0, x1, x2, t], x2, x2)/3)*u2[x0, x1, x2, t]*Re**(-1) - Derivative((p[x0, x1, x2, t] + rhoE[x0, x1, x2, t])*u0[x0, x1, x2, t], x0) - Derivative((p[x0, x1, x2, t] + rhoE[x0, x1, x2, t])*u1[x0, x1, x2, t], x1) - Derivative((p[x0, x1, x2, t] + rhoE[x0, x1, x2, t])*u2[x0, x1, x2, t], x2) - Derivative(q0[x0, x1, x2, t], x0) - Derivative(q1[x0, x1, x2, t], x1) - Derivative(q2[x0, x1, x2, t], x2)"


def test_maximum_derivative_order(mass, momentum, energy):
    """ Check that the maximum Derivative order is second-order. """
    assert maximum_derivative_order(flatten([mass.expanded, momentum.expanded, energy.expanded])) == 2


def test_remove_repeated_index():
    """ Check that repeated indices are removed. """
    indices = [Idx("i"), Idx("j"), Idx("j")]
    assert remove_repeated_index(indices) == [Idx("i")]


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))

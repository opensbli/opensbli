#!/usr/bin/env python

import os
import pytest

# OpenSBLI functions
from opensbli.grid import *
from opensbli.problem import *
from opensbli.utils import *

@pytest.fixture
def grid():
    return Grid(ndim=2)


@pytest.fixture
def coordinate_symbol():
    return "x"


@pytest.fixture
def mass(coordinate_symbol):
    return Equation("Eq(Der(rho,t), -Conservative(rhou_j,x_j))", 2, coordinate_symbol, substitutions=[], constants=[])


@pytest.fixture
def navier_stokes_problem(coordinate_symbol):
    # Number of dimensions
    ndim = 2

    # Equations
    mass = "Eq(Der(rho,t),- Conservative(rhou_j,x_j))"
    momentum = "Eq(Der(rhou_i,t) ,-Conservative(rhou_i*u_j + p* KD(_i,_j),x_j) + Der(tau_i_j,x_j) )"
    energy = "Eq(Der(rhoE,t),- Conservative((p+rhoE)*u_j,x_j) +Der(q_j,x_j) + Der(u_i*tau_i_j ,x_j) )"
    equations = [mass, momentum, energy]

    # Substitutions
    stress_tensor = "Eq(tau_i_j, (mu)*(Der(u_i,x_j)+ Conservative(u_j,x_i)- (2/3)* KD(_i,_j)* Der(u_k,x_k)))"
    heat_flux = "Eq(q_j,  (mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_j))"
    substitutions = [stress_tensor, heat_flux]

    # Define all the constants in the equations
    constants = ["Re", "Pr","gama","mu", "Minf", "C23", "c_j"]

    # Metrics
    metrics = [False, False]

    # Formulas for the variables used in the equations
    velocity = "Eq(u_i, rhou_i/rho)"
    pressure = "Eq(p, (gama-1)*(rhoE - (1/(2))*(u_j*u_j)))"
    temperature = "Eq(T, p*gama*Minf*Minf/(rho))"
    viscosity = "Eq(mu, T**(2/3))"
    formulas = [velocity, pressure, temperature, viscosity]

    problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
    return problem


def test_get_indexed_variables(mass, grid):
    """ Ensure that all of the Indexed grid variables are returned. """

    variables, count = get_indexed_variables([mass.expanded])
    
    rho = IndexedBase("rho")
    rhou0 = IndexedBase("rhou0")
    rhou1 = IndexedBase("rhou1")
    x0 = EinsteinTerm("x0")
    x1 = EinsteinTerm("x1")
    t = EinsteinTerm("t")

    assert variables == [rho[x0, x1, t], rhou0[x0, x1, t], rhou1[x0, x1, t]]

    return 


def test_get_derivatives(navier_stokes_problem, grid):
    """ Ensure that spatial and temporal Derivative objects are identified and returned correctly. """

    expanded_equations = navier_stokes_problem.get_expanded(navier_stokes_problem.equations)
    expanded_formulas = navier_stokes_problem.get_expanded(navier_stokes_problem.formulas)
    
    spatial_derivatives, temporal_derivatives = get_derivatives(flatten(expanded_equations))
    
    assert str(spatial_derivatives) == "[Derivative(rhou0[x0, x1, t], x0), Derivative(rhou1[x0, x1, t], x1), Derivative(rhou0[x0, x1, t]*u1[x0, x1, t], x1), Derivative(p[x0, x1, t] + rhou0[x0, x1, t]*u0[x0, x1, t], x0), Derivative(u1[x0, x1, t], x0, x1), Derivative(u0[x0, x1, t], x0, x0), Derivative(u0[x0, x1, t], x1, x1), Derivative(rhou1[x0, x1, t]*u0[x0, x1, t], x0), Derivative(p[x0, x1, t] + rhou1[x0, x1, t]*u1[x0, x1, t], x1), Derivative(u0[x0, x1, t], x0, x1), Derivative(u1[x0, x1, t], x1, x1), Derivative(u1[x0, x1, t], x0, x0), Derivative((p[x0, x1, t] + rhoE[x0, x1, t])*u0[x0, x1, t], x0), Derivative((p[x0, x1, t] + rhoE[x0, x1, t])*u1[x0, x1, t], x1), Derivative(u0[x0, x1, t], x0), Derivative(u1[x0, x1, t], x1), Derivative(u0[x0, x1, t], x1), Derivative(u1[x0, x1, t], x0), Derivative(T[x0, x1, t], x0, x0), Derivative(T[x0, x1, t], x1, x1)]"
    
    assert str(temporal_derivatives) == "[Derivative(rho[x0, x1, t], t), Derivative(rhou0[x0, x1, t], t), Derivative(rhou1[x0, x1, t], t), Derivative(rhoE[x0, x1, t], t)]"
    
    return 

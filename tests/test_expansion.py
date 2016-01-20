#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.expansion import *

def test_read_input():
    """ Ensure that the equation inputs are read correctly. """
    
    i = read_input(os.path.dirname(os.path.abspath(__file__)) + "/data/equations")
    
    # Equations
    assert len(i.equations) == 3
    assert i.equations == ["Eq(Der(rho,t),- conser(rhou_j,x_j))", "Eq(Der(rhou_i,t) ,-conser(rhou_i*u_j + p* KroneckerDelta(_i,_j),x_j) + Der(tau_i_j,x_j) )", "Eq(Der(rhoE,t),-conser((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )"]
    
    # Substitutions
    assert len(i.substitutions) == 2
    assert i.substitutions == ["Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j) + Der(u_j,x_i) - (2/3) * KroneckerDelta(_i,_j)*Der(u_k,x_k)))", "Eq(q_i, -(mu/((gama-1)*Minf*Minf*Pr*Re))*Der(T,x_i))"]

    # Dimension
    assert i.ndim == 3
    
    # Constants
    assert len(i.constants) == 6
    assert i.constants == ['Re', 'Pr', 'mu', 'gama', 'Minf', 'C23']

    # Coordinate symbol
    assert i.coordinate_symbol == ["x"]
    
    # Metrics
    assert i.metrics == ["False", "False"]
    
    # Formulas
    assert len(i.formulas) == 4
    assert i.formulas == ['Eq(u_i, rhou_i/rho)', 'Eq(p, (gama-1)*(rhoE - (1/(2*rho))*(rhou_j**2)))', 'Eq(T, p*gama*Minf*Minf/(rho))', 'Eq(mu, T**(2/3))']
    
    return

if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))

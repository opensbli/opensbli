#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.expansion import *

def test_read_input():
    """ Ensure that the equation inputs are read correctly. """
    
    i = read_input(os.path.dirname(os.path.abspath(__file__)) + "/data/equations")
    
    # Check equations
    assert len(i.equations) == 3
    assert i.equations == ["Eq(Der(rho,t),- conser(rhou_j,x_j))", "Eq(Der(rhou_i,t) ,-conser(rhou_i*u_j + p* KroneckerDelta(_i,_j),x_j) + Der(tau_i_j,x_j) )", "Eq(Der(rhoE,t),-conser((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )"]

    return

if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))

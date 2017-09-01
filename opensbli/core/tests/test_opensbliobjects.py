import os
import pytest
from sympy import pprint, Idx, srepr
from opensbli.core.opensbliobjects import *

@pytest.fixture
def scalar():
	return EinsteinTerm('u')

@pytest.fixture
def vector():
	return EinsteinTerm('u_i')

@pytest.fixture
def tensor():
	return EinsteinTerm('u_i_j')

@pytest.fixture
def spatial_coordinate():
	return CoordinateObject('x_i')

@pytest.fixture
def time_coordinate():
	return CoordinateObject('t', **{'time': True})

@pytest.fixture
def constant():
	return ConstantObject('a')

def test_ConstantObject(constant):
	# Check constant has no indices
	assert len(constant.indices) == 0
	# Check it is constant
	assert constant.is_constant == True
	# Check it is not a coordinate
	assert constant.is_coordinate == False
	# Check name
	assert constant.name == 'a'
	return

def test_CoordinateObject(spatial_coordinate, time_coordinate):
	# Check these objects are coordinates
	assert spatial_coordinate.is_coordinate == True
	assert time_coordinate.is_coordinate == True
	# Check indices
	assert len(spatial_coordinate.indices) == 1
	assert len(time_coordinate.indices) == 0
	# Check if time coordinate attribute
	assert spatial_coordinate.timecoordinate == False
	assert time_coordinate.timecoordinate == True
	return

def test_EinsteinTerm(scalar, vector, tensor):
	# Check properties of the EinsteinTerm
	assert scalar.name == 'u'
	assert vector.name == 'u_i'
	assert tensor.name == 'u_i_j'
	# Check indices
	assert len(scalar.get_indices()) == 0
	assert len(vector.get_indices()) == 1
	assert len(tensor.get_indices()) == 2
	# Check base names
	assert scalar.get_base() == 'u'
	assert vector.get_base() == 'u'
	assert tensor.get_base() == 'u'
	# Check indices are of the correct type
	assert type(scalar.get_indices()) == list
	assert type(vector.get_indices()[0]) == Idx
	assert type(tensor.get_indices()[0]) == Idx
	assert type(tensor.get_indices()[1]) == Idx
	# Check string represenation of the structure
	assert srepr(scalar.structure()) == "EinsteinTerm('u')"
	assert srepr(vector.structure()) == "Indexed(IndexedBase(Symbol('u')), Idx(Symbol('i', integer=True)))"
	assert srepr(tensor.structure()) == "Indexed(IndexedBase(Symbol('u')), Idx(Symbol('i', integer=True)), Idx(Symbol('j', integer=True)))"
	# Check indices are applied correctly
	assert srepr(vector.apply_index(vector.indices[0], 3)) == "EinsteinTerm('u3')"
	assert srepr(tensor.apply_index(tensor.indices[0], 4)) == "EinsteinTerm('u4_j')"
	# Check application of one index hasn't changed the other
	a = tensor.apply_index(tensor.indices[1], 4)
	assert a.indices[0] == tensor.indices[0]
	a = tensor.apply_index(tensor.indices[0], 5)
	assert a.indices[0] == tensor.indices[1]
	assert srepr(tensor.apply_index(tensor.indices[1], 5)) == "EinsteinTerm('u_i5')"
	# Check multi indices are applied correctly, ## check for python3 zip 
	assert srepr(tensor.apply_multiple_indices(tensor.indices, (zip(tensor.indices, [5,6])))) == "EinsteinTerm('u56')"
	return







if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
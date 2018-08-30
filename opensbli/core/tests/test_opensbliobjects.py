
#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (c) see License file

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

import os
import pytest
from sympy import Idx, srepr, Indexed
from opensbli.core.opensbliobjects import CoordinateObject, ConstantObject, Constant, EinsteinTerm,\
    ConstantIndexed, MetricObject, DataObject, DataSetBase, DataSet, IndexedBase
from opensbli.core.block import SimulationBlock


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


@pytest.fixture
def constant_indexed():
    return ConstantIndexed('a', Idx('i', 3))


@pytest.fixture
def metric():
    return MetricObject('xi')


@pytest.fixture
def data_obj():
    return DataObject('p')


@pytest.fixture
def local_block():
    return SimulationBlock(3, block_number=0)


@pytest.fixture
def datasetbase():
    return DataSetBase('rho')


def test_DataSet(local_block, datasetbase):
    dataset = datasetbase[[0, 0, 0]]
    # Check grid indices
    assert dataset.get_grid_indices == [0, 0, 0]
    # Check type
    assert isinstance(dataset, Indexed) is True
    return


def test_DataSetBase(local_block, datasetbase):
    # Check type
    assert isinstance(datasetbase, IndexedBase)
    # Check DataSetBase offset
    assert datasetbase._offset == 0
    # Check blockname and number
    assert datasetbase.blockname == local_block.blockname
    assert datasetbase.blocknumber == local_block.blocknumber
    # Check block independent name
    assert datasetbase.noblockname == EinsteinTerm('rho')
    # Check shape, for ndim = 3
    assert len(datasetbase._shape) == 4
    # Check simplelabel, SymPy string and hashable content
    assert datasetbase.simplelabel() == 'rho'
    # Should this have the p input argument?
    assert datasetbase._sympystr(None) == 'rho_B0'
    assert datasetbase._hashable_content() == 'rhoopensbliblock00'
    # Check default indices are zeros
    assert datasetbase.location() == [0, 0, 0]
    # Check DataSetBase converts to a DataSet
    dset = datasetbase[[1, 2, 3]]
    assert isinstance(dset, DataSet) is True
    return


def test_DataObject(data_obj):
    # Check DataObject is not a coordinate or constant
    assert data_obj.is_constant is False
    assert data_obj.is_coordinate is False
    # Check free_symbols
    assert DataObject('p') in data_obj.free_symbols
    assert type(data_obj.free_symbols) is set
    # Check DataObject is an EinsteinTerm
    assert isinstance(data_obj, EinsteinTerm) is True
    return


def test_MetricObject(metric):
    # Check metric is not a time coordinate
    assert metric.timecoordinate is False
    return


def test_ConstantIndexed(constant_indexed):
    # Check constant has an index
    assert len(constant_indexed.indices) is not 0
    # Check upper and lower bounds of the index
    assert constant_indexed.indices[0].lower == 0
    assert constant_indexed.indices[0].upper == 2
    # Check it is a constant
    assert constant_indexed.is_constant is True
    # Check type of index is Idx
    assert isinstance(constant_indexed.location[0], Idx) is True
    # Check object is Indexed and Constant
    assert isinstance(constant_indexed, Indexed) is True
    assert isinstance(constant_indexed, Constant) is True
    return


def test_ConstantObject(constant):
    # Check constant has no indices
    assert len(constant.indices) is 0
    # Check it is constant
    assert constant.is_constant is True
    # Check it is not a coordinate
    assert constant.is_coordinate is False
    # Check name and isinstance Constant class
    assert constant.name == 'a'
    assert isinstance(constant, Constant) is True
    return


def test_CoordinateObject(spatial_coordinate, time_coordinate):
    # Check these objects are coordinates
    assert spatial_coordinate.is_coordinate is True
    assert time_coordinate.is_coordinate is True
    # Check indices
    assert len(spatial_coordinate.indices) is 1
    assert len(time_coordinate.indices) is 0
    # Check if time coordinate attribute
    assert spatial_coordinate.timecoordinate is False
    assert time_coordinate.timecoordinate is True
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
    # No indices for scalar objects
    assert len(scalar.get_indices()) == 0
    ## Check all the indices of vector and tensor are type Idx
    assert [isinstance(ind, Idx) for ind in vector.get_indices()] == [True]
    assert [isinstance(ind, Idx) for ind in tensor.get_indices()] == [True]*2
    # Check string represenation of the structure used for Einstein expansion
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
    assert srepr(tensor.apply_multiple_indices(tensor.indices, (zip(tensor.indices, [5, 6])))) == "EinsteinTerm('u56')"
    return


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))

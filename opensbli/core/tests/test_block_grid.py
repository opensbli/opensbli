from opensbli.core.block import SimulationBlock
from sympy import srepr, pprint, Idx
from opensbli.core.opensbliobjects import ConstantObject
import pytest

@pytest.fixture
def block0_2d():
    return SimulationBlock(2)

def test_block_2d(block0_2d):
    assert block0_2d.ndim == 2
    assert len(block0_2d.Idxed_shape) == 2
    assert len(block0_2d.grid_indexes) == 2
    assert len(block0_2d.deltas) == 2
    assert len(block0_2d.shape) == 2
    assert [isinstance(t, Idx) for t in block0_2d.Idxed_shape] == [True, True]
    assert [t.upper for t in block0_2d.Idxed_shape] == block0_2d.shape
    assert [[t.lower, t.upper] for t in block0_2d.Idxed_shape] == block0_2d.ranges
    assert block0_2d.work_index == 0
    assert [isinstance(t, ConstantObject) for t in block0_2d.shape] == [True]*block0_2d.ndim
    assert block0_2d.stored_index == 0
    assert block0_2d.blocknumber == 0
    assert block0_2d.group_derivatives == False
    assert block0_2d.kernel_counter == 0
    return

def test_workaray_block0_2d(block0_2d):
    
    assert srepr(block0_2d.work_array()) == """DataSet(DataSetBase(Symbol('wk0'), Tuple(ConstantObject('block0np0', integer=True), ConstantObject('block0np1', integer=True), Idx(Symbol('opensbliblock00', integer=True)))), Integer(0), Integer(0), Idx(Symbol('opensbliblock00', integer=True)))"""
    block0_2d.increase_work_index  # Increases the work array index
    assert block0_2d.work_index == 1

    assert srepr(block0_2d.work_array()) == """DataSet(DataSetBase(Symbol('wk1'), Tuple(ConstantObject('block0np0', integer=True), ConstantObject('block0np1', integer=True), Idx(Symbol('opensbliblock00', integer=True)))), Integer(0), Integer(0), Idx(Symbol('opensbliblock00', integer=True)))"""
    block0_2d.reset_work_index  # Resets the work index

    assert srepr(block0_2d.work_array()) == """DataSet(DataSetBase(Symbol('wk0'), Tuple(ConstantObject('block0np0', integer=True), ConstantObject('block0np1', integer=True), Idx(Symbol('opensbliblock00', integer=True)))), Integer(0), Integer(0), Idx(Symbol('opensbliblock00', integer=True)))"""

    block0_2d.increase_work_index  # Increase the work array index to 1
    block0_2d.store_work_index  # Store the work array index (1)

    assert srepr(block0_2d.work_array()) == """DataSet(DataSetBase(Symbol('wk1'), Tuple(ConstantObject('block0np0', integer=True), ConstantObject('block0np1', integer=True), Idx(Symbol('opensbliblock00', integer=True)))), Integer(0), Integer(0), Idx(Symbol('opensbliblock00', integer=True)))"""
    block0_2d.increase_work_index  # Increase the index to 2

    assert srepr(block0_2d.work_array()) == """DataSet(DataSetBase(Symbol('wk2'), Tuple(ConstantObject('block0np0', integer=True), ConstantObject('block0np1', integer=True), Idx(Symbol('opensbliblock00', integer=True)))), Integer(0), Integer(0), Idx(Symbol('opensbliblock00', integer=True)))"""

    block0_2d.increase_work_index  # Increase once more to 3
    assert block0_2d.work_index == 3

    block0_2d.reset_work_to_stored  # resets the work array index to 1, previously stored
    assert srepr(block0_2d.work_array()) == """DataSet(DataSetBase(Symbol('wk1'), Tuple(ConstantObject('block0np0', integer=True), ConstantObject('block0np1', integer=True), Idx(Symbol('opensbliblock00', integer=True)))), Integer(0), Integer(0), Idx(Symbol('opensbliblock00', integer=True)))"""
    
    block0_2d.reset_work_index
    
    assert block0_2d.work_index == 0
    return

def test_kernel_counter_block(block0_2d):
    block0_2d.increase_kernel_counter
    assert block0_2d.kernel_counter == 1    
    block0_2d.increase_kernel_counter
    assert block0_2d.kernel_counter == 2
    block0_2d.reset_kernel_counter
    assert block0_2d.kernel_counter == 0
    return

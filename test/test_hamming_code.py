import pytest
from ucc_ft.hamming import HammingCode


def test_z_stabilizer_idx():
    # Test for r = 4, i.e. [15, 7, 3]
    hc = HammingCode(4)
    assert hc.z_stabilizer_idx(0) == list(range(7, 15))
    assert hc.z_stabilizer_idx(1) == list(range(3, 7)) + list(range(11, 15))
    assert hc.z_stabilizer_idx(2) == [1, 2, 5, 6, 9, 10, 13, 14]
    assert hc.z_stabilizer_idx(3) == list(range(0, 15, 2))
    with pytest.raises(AssertionError):
        hc.z_stabilizer_idx(5)


def test_x_stabilizer_idx():
    # Test for r = 4, i.e. [15, 7, 3]
    hc = HammingCode(4)
    assert hc.z_stabilizer_idx(0) == list(range(7, 15))
    assert hc.z_stabilizer_idx(1) == list(range(3, 7)) + list(range(11, 15))
    assert hc.z_stabilizer_idx(2) == [1, 2, 5, 6, 9, 10, 13, 14]
    assert hc.z_stabilizer_idx(3) == list(range(0, 15, 2))
    with pytest.raises(AssertionError):
        hc.z_stabilizer_idx(5)


def test_logical_x_idx():
    logical_x_idx = [
        [2, 4, 5],
        [2, 8, 9],
        [4, 6, 8, 10],
        [4, 8, 11],
        [2, 6, 8, 12],
        [6, 8, 13],
        [2, 4, 8, 14],
    ]
    assert HammingCode(4).logical_x_idx() == logical_x_idx


def test_logical_z_idx():
    logical_z_idx = [
        [1, 3, 5],
        [1, 7, 9],
        [0, 1, 7, 10],
        [3, 7, 11],
        [0, 3, 7, 12],
        [1, 3, 7, 13],
        [0, 1, 3, 7, 14],
    ]
    assert HammingCode(4).logical_z_idx() == logical_z_idx


def test_physical_z_idx():
    hc = HammingCode(4)
    assert sorted(hc.physical_z_idx()) == list(range(0, 15))

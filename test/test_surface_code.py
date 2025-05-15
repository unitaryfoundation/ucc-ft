import pytest
from ucc_ft.surface_code import RotatedSurfaceCode


def test_z_stabilizer():
    # Test for d=3
    sc = RotatedSurfaceCode(3)
    assert sc.z_stabilizer_idx(0) == [0, 1]
    assert sc.z_stabilizer_idx(1) == [1, 2, 4, 5]
    assert sc.z_stabilizer_idx(2) == [3, 4, 6, 7]
    assert sc.z_stabilizer_idx(3) == [7, 8]
    with pytest.raises(AssertionError):
        sc.z_stabilizer_idx(5)


def test_rotate():
    # Test for d=3
    orig = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    rotated = [[6, 3, 0], [7, 4, 1], [8, 5, 2]]
    sc = RotatedSurfaceCode(3)
    for i in range(3):
        for j in range(3):
            assert sc.rotate_idx(orig[i][j]) == rotated[i][j]


def test_x_stabilizer():
    # Test for d=3
    sc = RotatedSurfaceCode(3)
    assert sorted(sc.x_stabilizer_idx(0)) == [3, 6]
    assert sorted(sc.x_stabilizer_idx(1)) == [0, 1, 3, 4]
    assert sorted(sc.x_stabilizer_idx(2)) == [4, 5, 7, 8]
    assert sorted(sc.x_stabilizer_idx(3)) == [2, 5]
    with pytest.raises(AssertionError):
        sc.x_stabilizer_idx(5)


def test_logical_x_idx():
    # Test for d=3
    assert RotatedSurfaceCode(3).logical_x_idx() == [3, 4, 5]
    assert RotatedSurfaceCode(5).logical_x_idx() == [10, 11, 12, 13, 14]


def test_logical_z_idx():
    assert sorted(RotatedSurfaceCode(3).logical_z_idx()) == [1, 4, 7]


def test_physical_z_idx():
    # Test for d=3
    sc = RotatedSurfaceCode(3)
    assert sorted(sc.physical_z_idx()) == [0, 1, 2, 3, 4, 5, 6, 7, 8]

from ucc_ft.iceberg_code import IcebergCode


def test_z_stabilizer():
    # Test for m = 2, i.e. [4, 2, 2]
    ic = IcebergCode(2)
    assert ic.z_stabilizer_idx() == [0, 1, 2, 3]
 

def test_x_stabilizer():
    # Test for m = 2, i.e. [4, 2, 2]
    ic = IcebergCode(2)
    assert ic.x_stabilizer_idx() == [0, 1, 2, 3]
    

def test_logical_x_idx():
    assert IcebergCode(2).logical_x_idx() == [[1, 2], [1, 3]]
    assert IcebergCode(3).logical_x_idx() == [[1, 2], [1, 3], [1, 4], [1, 5]]


def test_logical_z_idx():
    assert IcebergCode(2).logical_z_idx() == [[0, 2], [0, 3]]
    assert IcebergCode(3).logical_z_idx() == [[0, 2], [0, 3], [0, 4], [0, 5]]


def test_physical_z_idx():
    ic = IcebergCode(2)
    assert sorted(ic.physical_z_idx()) == [0, 1, 2, 3]

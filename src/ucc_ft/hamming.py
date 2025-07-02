import numpy as np
from qldpc.objects import Pauli
from qldpc.codes import CSSCode
from stim import PauliString
from typing import List


class HammingCode(CSSCode):
    """Quantum error correcting code: [2^r - 1, 2^r - 2r - 1, 3]."""

    def __init__(self, size: int):
        self.size = size
        # Construct the r x 2^r - 1 parity check matrix
        # Each column is a unique nonzero binary vector of length r
        cols = [[int(x) for x in format(i, f"0{size}b")] for i in range(1, 2**size)]
        self.H = np.array(cols).T  # Shape: r x n
        super().__init__(code_x=self.H, code_z=self.H, field=2, is_subsystem_code=False)

    def z_stabilizer_idx(self, idx):
        """
        Return the qubit indices for the `idx`-th Z stabilizer for the Hamming code.
        """
        assert idx < self.size, "Index out of bounds"
        zs = self.get_stabilizer_ops(pauli=Pauli.Z)
        res = [i for i, z in enumerate(zs[idx]) if z]
        return res

    def x_stabilizer_idx(self, idx: int):
        """
        Return the qubit indices for the `idx`-th X stabilizer for the Hamming code.
        """
        assert idx < self.size, "Index out of bounds"
        xs = self.get_stabilizer_ops(pauli=Pauli.X)
        res = [i for i, x in enumerate(xs[idx]) if x]
        return res

    def stabilizers(self) -> List[PauliString]:
        """Return list of PauliString for the stabilizers of the
        iceberg code."""
        res = []
        for i in range(self.H.shape[0]):
            res.append(
                PauliString("*".join([f"Z{j}" for j in self.z_stabilizer_idx(i)])),
            )
            res.append(
                PauliString("*".join([f"X{j}" for j in self.x_stabilizer_idx(i)])),
            )
        return res

    def logical_x_idx(self) -> List[List[int]]:
        """Return the logical X operator for the iceberg code"""
        x_ops = self.get_logical_ops(pauli=Pauli.X)
        return [[i for i, x in enumerate(op) if x] for op in x_ops]

    def logical_z_idx(self) -> List[List[int]]:
        """Return the logical Z operator for the iceberg code."""
        z_ops = self.get_logical_ops(pauli=Pauli.Z)
        return [[i for i, z in enumerate(op) if z] for op in z_ops]

    def logical_x(self) -> List[PauliString]:
        """Return the logical X operator for the iceberg code."""
        res = []
        for x_ind in self.logical_x_idx():
            res.append(PauliString("*".join([f"X{j}" for j in x_ind])))
        return res

    def logical_z(self) -> List[PauliString]:
        """Return the logical Z operator for the iceberg code."""
        res = []
        for z_ind in self.logical_z_idx():
            res.append(PauliString("*".join([f"Z{j}" for j in z_ind])))
        return res

    def logical_prep_stabilizer(self) -> PauliString:
        """The prepared state is |0>_L, the +1 eigenstate of the logical Z operator."""
        return self.logical_z()

    def physical_z_idx(self) -> List[int]:
        """Return the physical Z operator for the iceberg code."""
        return [i for i in range(self.num_qubits)]

    def physical_z_stabilizers(self) -> List[PauliString]:
        """Return the physical Z operator for the iceberg code."""
        return [PauliString(f"Z{j}") for j in self.physical_z_idx()]

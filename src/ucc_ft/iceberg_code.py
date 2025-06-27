from qldpc.objects import Pauli
from qldpc.codes import IcebergCode
from stim import PauliString
from typing import List


class IcebergCode(IcebergCode):
    """Quantum error detecting code: [2m, 2m-2, 2].

    The m = 3 IcebergCode is the [6, 4, 2] code that is used to construct concatenated
    many-hypercube codes.

    References:
    - https://errorcorrectionzoo.org/c/iceberg
    - https://errorcorrectionzoo.org/c/stab_6_4_2
    - https://arxiv.org/abs/2403.16054
    """

    super.__init__()

    def z_stabilizer_idx(self):
        """
        Return the qubit indices for the Z stabilizer of the iceberg code.
        """
        zs = self.get_stabilizer_ops(pauli=Pauli.Z)
        res = [i for i, z in enumerate(zs[0]) if z]
        return res

    def x_stabilizer_idx(self):
        """
        Return the qubit indices for the idx-th Z stabilizer for rotated surface code
        of distance d.
        """
        # Rotate the Z stabilizer indices by 90 degrees clockwise
        xs = self.get_stabilizer_ops(pauli=Pauli.X)
        res = [i for i, x in enumerate(xs[0]) if x]
        return res

    def stabilizers(self) -> List[PauliString]:
        """Return list of PauliString for the stabilizers of the
        iceberg code."""
        res = [
            PauliString("*".join([f"Z{j}" for j in self.z_stabilizer_idx()])),
            PauliString("*".join([f"X{j}" for j in self.x_stabilizer_idx()])),
        ]
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

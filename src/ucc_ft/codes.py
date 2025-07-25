from stim import PauliString
from typing import List

# For the rotated surface code, we follow the convention in the provided sample code.
# The qubits are laid out on a grid, where qubits are numbered in row major order,
# and start with $X$ stabilizer in the top-left plaquette. The layout for
# the $d=3$ code is below:
# ```
#          Z
#     q0 •───• q1 ──• q2
#        │ X │  Z   |    X
#     q3 •───• q4 ──• q5
#  X     │ Z │  X   |
#     q6 •───• q7 ──• q8
#                 Z
# ```
#
# Note that rotating the above layout by 90 degrees swaps X and Z stabilizers


class RotatedSurfaceCode:
    def __init__(self, d: int):
        if d % 2 == 0:
            raise ValueError("Distance d must be odd for a rotated surface code.")
        self.d = d

    @property
    def num_qubits(self) -> int:
        """Return the number of qubits in the rotated surface code of distance d."""
        return self.d * self.d

    def z_stabilizer_idx(self, idx: int):
        """
        Return the qubit indices for the idx-th Z stabilizer for rotated surface code
        of distance d.
        """
        assert idx < (self.d * self.d - 1) // 2, "Index out of bounds"

        if idx < (self.d - 1) // 2:
            return [idx * 2, idx * 2 + 1]
        elif idx >= self.d * (self.d - 1) // 2:
            return [2 * idx + 1, 2 * (idx + 1)]
        else:
            i = idx // ((self.d - 1) // 2)
            j = ((idx % ((self.d - 1) // 2)) * 2) + 1 + (i % 2)
            return [
                (i - 1) * self.d + j - 1,
                (i - 1) * self.d + j,
                i * self.d + j - 1,
                i * self.d + j,
            ]

    def rotate_idx(self, idx: int):
        """Given a qubit index in the row-major order, return the corresponding index
        if it were rotated by 90 degrees clockwise.
        """
        i = idx // self.d + 1
        j = idx % self.d + 1
        return (self.d - j) * self.d + i - 1

    def x_stabilizer_idx(self, idx: int):
        """
        Return the qubit indices for the idx-th Z stabilizer for rotated surface code
        of distance d.
        """
        # Rotate the Z stabilizer indices by 90 degrees clockwise
        zs = self.z_stabilizer_idx(idx)
        res = []
        for z_idx in zs:
            res.append(self.rotate_idx(z_idx))
        return res

    def stabilizers(self) -> List[PauliString]:
        """Return list of PauliString for the stabilizers of the
        rotated surface code of distance d."""
        res = []
        for i in range((self.d * self.d - 1) // 2):
            res.append(
                PauliString("*".join([f"Z{j}" for j in self.z_stabilizer_idx(i)]))
            )
            res.append(
                PauliString("*".join([f"X{j}" for j in self.x_stabilizer_idx(i)]))
            )
        return res

    def logical_x_idx(self) -> List[int]:
        """Logical X operator, go across the middle of the surface code to match paper"""
        offset = self.d * (self.d - 1) // 2
        return [offset + i for i in range(self.d)]

    def logical_z_idx(self) -> List[int]:
        """Return the logical Z operator for the rotated surface code of distance d."""
        return [self.rotate_idx(i) for i in self.logical_x_idx()]

    def logical_x(self) -> PauliString:
        """Return the logical X operator for the rotated surface code of distance d."""
        return PauliString("*".join([f"X{j}" for j in self.logical_x_idx()]))

    def logical_z(self) -> PauliString:
        """Return the logical Z operator for the rotated surface code of distance d."""
        return PauliString("*".join([f"Z{j}" for j in self.logical_z_idx()]))

    def logical_prep_stabilizer(self) -> PauliString:
        """The prepared state is |0>_L, the +1 eigenstate of the logical Z operator."""
        return self.logical_z()

    def physical_z_idx(self) -> List[int]:
        """Return the physical Z operator for the rotated surface code of distance d."""
        return [i for i in range(self.d * self.d)]

    def physical_z_stabilizers(self) -> List[PauliString]:
        """Return the physical Z operator for the rotated surface code of distance d."""
        return [PauliString(f"Z{j}") for j in self.physical_z_idx()]


class CatStateCode:
    """
    Although not an error-correcting code, there are a set of stabilizers and
    logical operators that are stabilized by the cat state. We can use this to
    verify the fault tolerance of the cat state preparation circuit.
    """

    def __init__(self, num_qubits: int, max_faults: int):
        self.num_qubits = num_qubits
        self.d = max_faults * 2 + 1

    def stabilizers(self):
        return [
            PauliString(f"Z{i}*Z{j}")
            for (i, j) in zip(range(self.num_qubits), range(1, self.num_qubits))
        ]

    def logical_prep_stabilizer(self):
        """The prepared state is |+>_L, the +1 eigenstate of the logical X operator."""
        return PauliString("X" * self.num_qubits)

    def physical_z_stabilizers(self):
        return [PauliString(f"Z{i}") for i in range(self.num_qubits)]

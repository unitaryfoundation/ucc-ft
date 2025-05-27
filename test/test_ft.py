import pytest
import numpy as np
from stim import PauliString, Tableau
from pathlib import Path  # Add this import

from ucc_ft.checker import to_julia_tableau_fmt
from ucc_ft.checker import (
    julia_source_to_qprog,
    ft_check_ideal,
    qasm_to_qprog,
    qasm_to_qprog_source,
)
from ucc_ft.surface_code import RotatedSurfaceCode


def test_to_julia_tableau_fmt():
    # 3-qubit bit flip example
    t = Tableau.from_stabilizers(
        [PauliString("ZZI"), PauliString("IZZ"), PauliString("XXX")]
    )
    res = to_julia_tableau_fmt(t)
    expected = np.array(
        [
            [
                False,
                False,
                False,
                True,
                True,
                False,
            ],
            [False, False, False, False, True, True],
            [True, True, True, False, False, False],
        ],
        dtype=bool,
    )
    assert np.array_equal(res, expected)


# class for "Cat State" code
class CatStateCode:
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


@pytest.mark.parametrize(
    "max_faults, expected_result",
    [
        (2, True),
        (3, False),
    ],
)
def test_check_ft_cat_state_with_different_faults(max_faults, expected_result):
    """Test the fault tolerance of the cat state circuit with varying max_faults."""

    num_qubits = 8
    code = CatStateCode(num_qubits, max_faults)

    circuit = """
    OPENQASM 3.0;
    include "stdgates.inc";

    const uint size = __NUM_QUBITS__;
    qubit[size] state;
    qubit ancilla;

    def cat_prep() {

        bit res = 1;
        while(res != 0) {
            reset state[0];
            res = 0;
            h state[0];

            // QASM ranges are inclusive for both start and end
            for int i in [1:(size-1)] {
                reset state[i];
                cx state[0], state[i];
            }

            // Parity check
            for int i in [1:(size-1)] {
                reset ancilla;
                cx state[i-1], ancilla;
                cx state[i], ancilla;
                bit tmp = measure ancilla;
                res = res | tmp;
            }
        }
    }
    """.replace("__NUM_QUBITS__", str(num_qubits))

    qprog_context = qasm_to_qprog(circuit)
    result = ft_check_ideal(
        code, qprog_context.get_qprog("cat_prep"), qprog_context, "prepare", NERRS=12
    )
    assert result == expected_result


def test_check_rotated_surface_prep():
    """Test the fault tolerance of the rotated surface code preparation circuit."""

    d = 3
    sc = RotatedSurfaceCode(d)
    julia_source_path = Path(__file__).parent / "rotated_surface_code.jl"
    julia_source = julia_source_path.read_text()

    qprog_context = julia_source_to_qprog(
        julia_source,
        [
            "_rotated_surface_prepare_0",
            "rotated_surface_prepare_0",
            "rotated_surface_z_m",
            "rotated_surface_x_m",
            "rotated_surface_lz_m",
            "_xadj",
            "_zadj",
            "prepare_cat",
            "prepare_cat_for_z",
            "generate_cat_verification",
            "mwpm_full",
        ],
    )
    result = ft_check_ideal(
        sc,
        qprog_context.get_qprog("rotated_surface_prepare_0", d),
        qprog_context,
        "prepare",
        NERRS=12,
    )
    assert result


def test_check_rotated_surface_prep_qasm():
    """Test the fault tolerance of the rotated surface code preparation circuit."""

    d = 3
    sc = RotatedSurfaceCode(d)
    qasm_source_path = Path(__file__).parent / "rotated_surface_code.qasm"
    qasm_source = qasm_source_path.read_text()
    qprog_src = qasm_to_qprog_source(qasm_source)

    julia_source_path = Path(__file__).parent / "rotated_surface_code.in_translation.jl"
    julia_source = julia_source_path.read_text()

    qprog_context = julia_source_to_qprog(
        julia_source + "\n\n" + qprog_src,
        [
            "prepare_state",
            "rotated_surface_z_m",
            "rotated_surface_x_m",
            "rotated_surface_lz_m",
            "_xadj",
            "_zadj",
            "prepare_cat_internal",
            "prepare_cat",
            "generate_cat_verification",
            "mwpm_full",
            "mwpm_full_x",
            "mwpm_full_z",
            "data_size",
            "cat_size",
            "state",
            "cat",
            "verify",
            "num_syndromes",
        ],
    )

    result = ft_check_ideal(
        sc,
        qprog_context.get_qprog("prepare_state"),
        qprog_context,
        "prepare",
        NERRS=12,
    )
    assert result


def test_check_rotated_surface_CNOT():
    """Test the fault tolerance of the rotated surface code gate circuit."""
    d = 3
    sc = RotatedSurfaceCode(d)

    circuit = """
    OPENQASM 3.0;
    include "stdgates.inc";

    const uint d = __d__;
    const uint data_size = d * d;
    qubit[data_size] state1;
    qubit[data_size] state2;

    def logical_CNOT() {
        // QASM ranges are inclusive for both start and end
        for int i in [0:(data_size-1)] {
            cx state1[i], state2[i];
        }
    }
    """.replace("__d__", str(d))
    qprog_context = qasm_to_qprog(circuit)

    result = ft_check_ideal(
        sc,
        qprog_context.get_qprog("logical_CNOT"),
        qprog_context,
        "gate",
        NERRS=12,
    )
    assert result


def test_check_rotated_surface_decoder():
    """Test the fault tolerance of the rotated surface code gate decoder+correction gadget."""
    d = 3
    sc = RotatedSurfaceCode(d)
    julia_source_path = Path(__file__).parent / "rotated_surface_code.jl"
    julia_source = julia_source_path.read_text()

    qprog_context = julia_source_to_qprog(
        julia_source,
        [
            "_rotated_surface_decoder",
            "rotated_surface_decoder",
            "rotated_surface_z_m",
            "rotated_surface_x_m",
            "_xadj",
            "_zadj",
            "prepare_cat",
            "generate_cat_verification",
            "mwpm",
            "mwpm2",
        ],
    )
    result = ft_check_ideal(
        sc,
        qprog_context.get_qprog("rotated_surface_decoder", d),
        qprog_context,
        "decoder",
        NERRS=12,
    )
    assert result


def test_check_rotated_surface_measure():
    """Test the fault tolerance of the rotated surface code gate measurement gadget."""
    d = 3
    sc = RotatedSurfaceCode(d)
    julia_source_path = Path(__file__).parent / "rotated_surface_code.jl"
    julia_source = julia_source_path.read_text()

    qprog_context = julia_source_to_qprog(
        julia_source,
        [
            "rotated_surface_measurement",
            "_rotated_surface_decoder",
            "rotated_surface_z_m",
            "rotated_surface_x_m",
            "rotated_surface_lz_m",
            "majority",
            "_xadj",
            "_zadj",
            "prepare_cat",
            "generate_cat_verification",
            "mwpm",
        ],
    )
    result = ft_check_ideal(
        sc,
        qprog_context.get_qprog("rotated_surface_measurement", d),
        qprog_context,
        "measurement",
        NERRS=12,
    )
    assert result

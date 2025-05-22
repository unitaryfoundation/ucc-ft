import pytest
import numpy as np
from stim import PauliString, Tableau
from pathlib import Path  # Add this import

from ucc_ft.checker import to_julia_tableau_fmt
from ucc_ft.checker import (
    julia_source_to_qprog,
    ft_check_ideal,
    qasm_to_qprog,
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

    # target_stabilizers = [
    #     PauliString(f"Z{i}*Z{j}")
    #     for (i, j) in zip(range(num_qubits), range(1, num_qubits))
    # ]
    # target_stabilizers.append(PauliString("X" * num_qubits))

    # input_stabilizers = [PauliString(f"Z{i}") for i in range(num_qubits)]

    circuit = """
    OPENQASM 3.0;
    include "stdgates.inc";

    const uint size = __NUM_QUBITS__;
    const uint ancilla_size = 1;
    def cat_prep(qubit[size] state, qubit[ancilla_size] ancilla) {

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
                reset ancilla[0];
                cx state[i-1], ancilla[0];
                cx state[i], ancilla[0];
                bit tmp = measure ancilla[0];
                res = res | tmp;
            }
        }
    }
    """.replace("__NUM_QUBITS__", str(num_qubits))

    result = ft_check_ideal(code, qasm_to_qprog(circuit), "prepare", NERRS=12)
    assert result == expected_result


def test_check_rotated_surface_prep():
    """Test the fault tolerance of the rotated surface code preparation circuit."""

    d = 3
    sc = RotatedSurfaceCode(d)
    julia_source_path = Path(__file__).parent / "rotated_surface_code.jl"
    julia_source = julia_source_path.read_text()

    qprog_and_context = julia_source_to_qprog(
        julia_source,
        "rotated_surface_prepare_0",
        [
            "_rotated_surface_prepare_0",
            "rotated_surface_prepare_0",
            "rotated_surface_z_m",
            "rotated_surface_x_m",
            "rotated_surface_lz_m",
            "_xadj",
            "_zadj",
            "prepare_cat",
            "generate_cat_verification",
            "mwpm_full",
        ],
    )
    qprog_and_context.qprog = qprog_and_context.qprog(d)
    result = ft_check_ideal(sc, qprog_and_context, "prepare", NERRS=12)
    assert result


def test_check_rotated_surface_CNOT():
    """Test the fault tolerance of the rotated surface code gate circuit."""
    d = 3
    sc = RotatedSurfaceCode(d)

    circuit = """
    OPENQASM 3.0;
    include "stdgates.inc";

    const uint size = __NUM_QUBITS__;

    def rotated_surface_CNOT(qubit[size] state1, qubit[size] state2) {

        // QASM ranges are inclusive for both start and end
        for int i in [0:(size-1)] {
            cx state1[i], state2[i];
        }

    }
    """.replace("__NUM_QUBITS__", str(sc.num_qubits))

    qprog_and_context = qasm_to_qprog(circuit)
    result = ft_check_ideal(sc, qprog_and_context, "gate", NERRS=12)
    assert result


def test_check_rotated_surface_decoder():
    """Test the fault tolerance of the rotated surface code gate decoder+correction gadget."""
    d = 3
    sc = RotatedSurfaceCode(d)
    julia_source_path = Path(__file__).parent / "rotated_surface_code.jl"
    julia_source = julia_source_path.read_text()

    qprog_and_context = julia_source_to_qprog(
        julia_source,
        "rotated_surface_decoder",
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
    qprog_and_context.qprog = qprog_and_context.qprog(d)
    result = ft_check_ideal(sc, qprog_and_context, "decoder", NERRS=12)
    assert result


def test_check_rotated_surface_measure():
    """Test the fault tolerance of the rotated surface code gate measurement gadget."""
    d = 3
    sc = RotatedSurfaceCode(d)
    julia_source_path = Path(__file__).parent / "rotated_surface_code.jl"
    julia_source = julia_source_path.read_text()

    qprog_and_context = julia_source_to_qprog(
        julia_source,
        "rotated_surface_measurement",
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
    qprog_and_context.qprog = qprog_and_context.qprog(d)
    result = ft_check_ideal(sc, qprog_and_context, "measurement", NERRS=12)
    assert result

    # TODO:
    #  Goal is ft_check(code: StabilizerCode, circuit: CircuitType, gadget_type: GadgetType)
    #           input_state = input_for(code, gadget_type) <--- can try?
    #           output_state = run_circuit(input_state) #<--messy for now, do by hand!
    #           # run the FT check
    #   Add other gadgets
    #    Stablizer code -> tableau (Where to add phases?)
    #    C
    # --- Custom QPROG
    #    1. Write QASM, but have external calls or predefined submodules for some elements?
    #         No multiqubit measur, but yes for cat state prep?
    #

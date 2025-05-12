import pytest
import numpy as np
from stim import PauliString, Tableau
from pathlib import Path  # Add this import

from ucc_ft.checker import to_julia_tableau_fmt
from ucc_ft.checker import ft_check, ft_check_from_qprog, julia_source_to_qprog
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

    target_stabilizers = [
        PauliString(f"Z{i}*Z{j}")
        for (i, j) in zip(range(num_qubits), range(1, num_qubits))
    ]
    target_stabilizers.append(PauliString("X" * num_qubits))

    input_stabilizers = [PauliString(f"Z{i}") for i in range(num_qubits)]

    circuit = """
    OPENQASM 3.0;
    include "stdgates.inc";

    const uint size = __NUM_QUBITS__;

    def cat_prep(qubit[size] state, qubit ancilla) {

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

    # Check if the circuit is fault tolerant
    d = max_faults * 2 + 1
    result = ft_check(input_stabilizers, target_stabilizers, circuit, d, NERRS=12)
    assert result == expected_result


def test_check_rotated_prep():
    """Test the fault tolerance of the rotated surface code preparation circuit."""

    d = 3
    sc = RotatedSurfaceCode(d)

    # The output of preparation is the logical-|0>, which requires
    # adding the logical Z operator to the list of stabilizers
    target = sc.stabilizers()
    target.append(sc.logical_z())

    # The initial input is all physical qubits in |0> state, which
    # requires adding the joint physical Z operator on all qubits
    # to the list of stabilizers
    init = sc.stabilizers()
    init.append(PauliString("*".join([f"Z{i}" for i in range(d * d)])))

    # Use pathlib to load the file from the resources directory
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
    result = ft_check_from_qprog(init, target, qprog_and_context, d, NERRS=12)
    assert result

    # Todo:
    # -- No custom QPROG (use as is)
    #    1. Generalize the check_FT call to take target and destination
    #               ( longer term comes from running the circuit?)
    #    2. Import the QPROG by hand
    #    3. Do the assert ... how to have a failure case too?
    # --- Custom QPROG
    #    1. Write QASM, but have external calls or predefined submodules for some elements?
    #         No multiqubit measur, but yes for cat state prep?
    #

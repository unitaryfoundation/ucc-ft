// Rotated surface code gadgets written in QASM
OPENQASM 3.0;
include "stdgates.inc";

const uint d = 3; // hard-coded for now
const uint size = d * d;

def logical_CNOT(qubit[size] state1, qubit[size] state2) {
    // QASM ranges are inclusive for both start and end
    for int i in [0:(size-1)] {
        cx state1[i], state2[i];
    }
}


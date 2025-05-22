// Rotated surface code gadgets written in QASM
OPENQASM 3.0;
include "stdgates.inc";

const uint d = 3; // hard-coded for now
const uint data_size = d * d;
const uint ancilla_size = d + 1; // should it be d -1 ??
const uint num_syndromes = (d*d-1)/2;


///////////////////////////////////////////////////////////////////
/// State Preparation
//  The "standard" state preparation where the state is initilialized,
//  and stabilizers are repeatedly measured until the state stabilizes.
//  Afterwards, the code state is corrected in to the logical-0 state based
//  on the MWPM matched errors


// Call into a classical MWPM function to determine if error occurred
extern mwpm_full_x(uint, bit[num_syndromes]) -> bit[data_size];
extern mwpm_full_z(uint, bit[num_syndromes], bit) -> bit[data_size];

def prepare_state(qubit[data_size] state, qubit[ancilla_size] ancilla) {

    for int i in [0:(data_size-1)] {
        reset state[i];
    }

    t = (d-1)/2 + 1;

    bit res = 1;

    bit[num_syndromes] s_x;
    bit[num_syndromes] s_z;
    bit s_lz = 0;

    while (res != 0) {
        res = 0;

        // Measure stabilizers
        for int round in [0:(t-1)] {
            for int j in [0:(num_syndromes-1)]{

                // HACK: +1 in second argument to match julia expected indexing!!
                bit m_x = rotated_surface_x_m(d, j+1);
                bit m_z = rotated_surface_z_m(d, j+1);

                // Check parity across rounds
                if(round > 0) {
                    res = res | (m_x ^ s_x[j]) | (m_z ^ s_z[j]);
                }
                // Update to latest measurement results
                s_x[j] = m_x;
                s_z[j] = m_z;
            }

            bit m_lz = rotated_surface_lz_m(d);
            if(round > 0) {
                res = res | (m_lz ^ s_lz);
            }
            s_lz = m_lz;
        }
        // res will now be 0 if all measurements were the same between rounds
    }

    // HACK!! mwpm_full expects to take Z3 context as first argument in Julia
    // For now, let that pass through ... fix this later
    bit[data_size] r_x = mwpm_full_x(ctx, d, s_x);
    bit[data_size] r_z = mwpm_full_z(ctx, d, s_z, s_lz);

    // Correct the MWPM marked errors
    for int i in [0:(data_size-1)] {
        if(r_x[i]) {
            z state[i];
        }
        if(r_z[i]) {
            x state[i];
        }
    }

    // Reset the ancilla
    for int i in [0:(ancilla_size-1)] {
        reset ancilla[i];
    }
}


/// Two-qubit gates
// Simple transversal application
def logical_CNOT(qubit[data_size] state1, qubit[data_size] state2) {
    // QASM ranges are inclusive for both start and end
    for int i in [0:(data_size-1)] {
        cx state1[i], state2[i];
    }
}


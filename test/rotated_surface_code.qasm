// Rotated surface code gadgets written in QASM
OPENQASM 3.0;
include "stdgates.inc";

const uint d = 3; // hard-coded for now
const uint data_size = d * d;
const uint cat_size = d+1;
const uint verify_size = 1;
const uint num_syndromes = (d*d-1)/2;

// TODO -- rewrite with qubits as d x d grid
qubit[data_size] state;
qubit[cat_size] cat;
qubit verify;

////////////////////////////////////////////////////////////////////////////////
// Conventions
// For the rotated surface code, we follow the convention in the provided sample code.
// The qubits are laid out on a grid, where qubits are numbered in row major order,
// and start with $X$ stabilizer in the top-left plaquette. The layout for
// the $d=3$ code is below:
// ```
//          Z
//     q0 •───• q1 ──• q2
//        │ X │  Z   |    X
//     q3 •───• q4 ──• q5
//  X     │ Z │  X   |
//     q6 •───• q7 ──• q8
//                 Z
// ```
//
// Note that rotating the above layout by 90 degrees swaps X and Z stabilizers

////////////////////////////////////////////////////////////////////////////////
//  State Preparation
//  The "standard" state preparation where the state is initilialized,
//  and stabilizers are repeatedly measured until the state stabilizes.
//  Afterwards, the code state is corrected in to the logical-0 state based
//  on the MWPM matched errors


// Call into a classical MWPM function to determine if error occurred
extern mwpm_full_x(uint, bit[num_syndromes]) -> bit[data_size];
extern mwpm_full_z(uint, bit[num_syndromes], bit) -> bit[data_size];

////////////////////////////////////////////////////////////////////////////////
//  Stabilizer measurements
//  This implementation follows the Shor method. For each stabilizer of
//  weight W, we prepare a cat state over W qubits, do a controlled-X (or -Z)
//  operation between corresponding data and cat state qubits, then measure the
//  ancilla qubits. This approach matches that in the reference Julia code, but
//  differs from more cannonical surface code measurement schecms that use say 1
//  measurement qubit per stabilizer (or some that re-use). Future versions will
//  consider proving fault tolerance for those variants.

// FT cat state preparation by verifying parity via extra qubit
//   Uses num_cat of cat_size qubits
def prepare_cat(uint num_cat) {
    bit res = 1;
    while(res != 0) {
        reset cat[0];
        res = 0;
        h cat[0];

        for int i in [1:(num_cat-1)] {
            reset cat[i];
            cx cat[0], cat[i];
        }

        for int i in [1:(num_cat-1)] {
            reset verify;
            cx cat[i-1], verify;
            cx cat[i], verify;
            bit tmp = measure verify;
            res = res | tmp;
        }
    }
}

// Measure the i-th X stabilizer
/*
def rotated_surface_x_m(uint idx) -> bit {

    uint num_cat = 2;
    prepare_cat_state(cat, verify);

    // boundaries
    for i in [0:X] {
        CNOT(cat[i], state[b[i]]);
    }
    // bulk

    // Measure ancilla in X basis to extract parity
    bit res = 0;
    for i in [0:(num_cat-1)] {
        h cat[i];
        bit tmp = measure cat[i];
        h cat[i];
        res = res ^ tmp;
    }
    return res;
}*/

// Measure the i-th Z stabilizer
def rotated_surface_z_m(uint idx) -> bit {
    uint num_cat = 2;

    if (idx < (d - 1) / 2) {
        num_cat = 2;
        prepare_cat(num_cat);

        cz cat[0], state[2 * idx];
        cz cat[1], state[2 * idx + 1];
    }
    if (idx >= d * (d-1)/ 2) {
        num_cat = 2;
        prepare_cat(num_cat);

        cz cat[0], state[2 * idx + 1 ];
        cz cat[1], state[2 * (idx + 1)];
    }
    if ((idx >= (d-1)/2) && (idx < d * (d-1)/2) ) {
        num_cat = 4;
        prepare_cat(num_cat);

        uint i = idx / ((d - 1) / 2);
        uint j = ((idx % ((d - 1) / 2)) * 2) + 1 + (i % 2);

        cz cat[0], state[(i - 1) * d + j - 1];
        cz cat[1], state[(i - 1) * d + j];
        cz cat[2], state[i * d + j - 1];
        cz cat[3], state[i * d + j];
    }

    // Measure ancilla in X basis to extract parity
    bit res = 0;
    for int c in [0:(num_cat-1)] {
        h cat[c];
        bit tmp = measure cat[c];
        h cat[c];
        res = res ^ tmp;
    }
    return res;
}

def prepare_state() {

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
                bit m_z = rotated_surface_z_m(j);

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
    for int i in [0:(cat_size-1)] {
        reset cat[i];
    }
    reset verify;
}


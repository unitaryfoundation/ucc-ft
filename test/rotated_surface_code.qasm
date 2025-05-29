// Rotated surface code gadgets written in QASM
OPENQASM 3.0;
include "stdgates.inc";

////////////////////////////////////////////////////////////////////////////////
// Conventions
// For the rotated surface code, we follow the convention in the Julia sample code.
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
// Code parameters
const uint d = 3;
const uint data_size = d * d;
const uint cat_size = d+1;
const uint verify_size = 1;
const uint num_syndromes = (d*d-1)/2;

////////////////////////////////////////////////////////////////////////////////
// Quantum state registers - QASM3 spec requires these be defined as globals
// TODO -- Rewrite with qubits as d x d grid versus flattened register
qubit[data_size] state; // the logical qubit
qubit[cat_size] cat;    // qubits prepared in cat state for syndrome measurement
qubit verify;           // used to FT verify cat state preparation

////////////////////////////////////////////////////////////////////////////////
//  State Preparation
//  The "standard" state preparation where stabilizers are repeatedly measured
//  until the state stabilizes. Afterwards, the code state is corrected in to
//  the logical-0 state based on the MWPM matched error.


////////////////////////////////////////////////////////////////////////////////
// External classical subroutines

// Calls into a classical MWPM function to determine if error occurred
// Takes in distance, syndrome measurement outcomes and returns whether to apply correction for that syndrome
extern mwpm_full_x(uint, bit[num_syndromes]) -> bit[data_size];
extern mwpm_full_z(uint, bit[num_syndromes], bit) -> bit[data_size];

// Given a row-major index into the surface grid, returns the corresponding
// row-major index if the grid were rotated 90 degrees
//   For d = 3 maps [[0,1,2],[3,4,5],[6,7,8]] to [[6,3,0],[7,4,1],[8,5,2]
extern rotate(uint, uint ) -> uint; //  d, index_in -> index_out

// For now, rely on this being extern to simplify translation to @qprog.
// That is, all other subroutines here are translated to julia @qprog, but this
// one shouldn't be since it just transform classical index to index value.
// Converting it to qprog adds unneeded symbolic execution (and doesn't quite work)

////////////////////////////////////////////////////////////////////////////////
//  Stabilizer measurements
//
//  This implementation follows the Shor method. For each stabilizer of
//  weight W, we prepare a cat state over W qubits, do a controlled-X (or -Z)
//  operation between corresponding data and cat state qubits, then measure the
//  parity of the cat qubits. This approach matches that in the reference Julia code, but
//  differs from more cannonical surface code measurement schemes that use 1
//  measurement qubit per stabilizer. Future versions will consider proving fault
//  tolerance for those variants.

// FT cat state preparation, using num_cat of the `cat` qubit register qubits
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

// Measure the i-th Z stabilizer
def rotated_surface_z_m(uint idx) -> bit {
    uint num_cat = 2;

    // Formulas here give the row-major index for that stabilizer

    // Top
    if (idx < (d - 1) / 2) {
        num_cat = 2;
        prepare_cat(num_cat);

        cz cat[0], state[2 * idx];
        cz cat[1], state[2 * idx + 1];
    }

    // Bottom
    if (idx >= d * (d-1)/ 2) {
        num_cat = 2;
        prepare_cat(num_cat);

        cz cat[0], state[2 * idx + 1 ];
        cz cat[1], state[2 * (idx + 1)];
    }

    // Inside
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

    // Measure cat state in X basis to extract parity
    bit res = 0;
    for int c in [0:(num_cat-1)] {
        h cat[c];
        bit tmp = measure cat[c];
        h cat[c];
        res = res ^ tmp;
    }
    return res;
}

// Measure the i-th X stabilizer
def rotated_surface_x_m(uint idx) -> bit {
    uint num_cat = 2;

    // Indices are `rotate` on the Z stabilizer indices above
    if (idx < (d - 1) / 2) {
        num_cat = 2;
        prepare_cat(num_cat);

        cx cat[0], state[rotate(d, 2 * idx)];
        cx cat[1], state[rotate(d, 2 * idx + 1)];
    }
    // Right-hand side
    if (idx >= d * (d-1)/ 2) {
        num_cat = 2;
        prepare_cat(num_cat);

        cx cat[0], state[rotate(d, 2 * idx + 1)];
        cx cat[1], state[rotate(d, 2 * (idx + 1))];
    }
    // Inside
    if ((idx >= (d-1)/2) && (idx < d * (d-1)/2) ) {
        num_cat = 4;
        prepare_cat(num_cat);

        uint i = idx / ((d - 1) / 2);
        uint j = ((idx % ((d - 1) / 2)) * 2) + 1 + (i % 2);

        cx cat[0], state[rotate(d, (i - 1) * d + j - 1)];
        cx cat[1], state[rotate(d, (i - 1) * d + j)];
        cx cat[2], state[rotate(d, i * d + j - 1)];
        cx cat[3], state[rotate(d, i * d + j)];
    }

    // Measure cat state in X basis to extract parity
    bit res = 0;
    for int c in [0:(num_cat-1)] {
        h cat[c];
        bit tmp = measure cat[c];
        h cat[c];
        res = res ^ tmp;
    }
    return res;
}



// Measure the Logical-Z operator
def rotated_surface_lz_m() -> bit {
    prepare_cat(d);

    // go down the middle set of qubits
    for int i in [0:(d-1)] {
        cz cat[i], state[(d*(2*i+1)-1)/2];
    }

    // Measure cat state in X basis to extract parity
    bit res = 0;
    for int c in [0:(d-1)] {
        h cat[c];
        bit tmp = measure cat[c];
        h cat[c];
        res = res ^ tmp;
    }
    return res;
}

// Prepare the logical-|0> state
def prepare_state() {

    for int i in [0:(data_size-1)] {
        reset state[i];
    }


    t = (d-1)/2 + 1; // number of rounds

    bit res = 1; // whether syndrome measurements have stabilized

    // syndrome outcomes at final round
    bit[num_syndromes] s_x;
    bit[num_syndromes] s_z;
    bit s_lz = 0;

    // repeat until success
    while (res != 0) {
        res = 0;


        for int round in [0:(t-1)] {

            // Measure X, Z stabilizers
            for int j in [0:(num_syndromes-1)]{


                bit m_x = rotated_surface_x_m(j);
                bit m_z = rotated_surface_z_m(j);

                // Check parity across rounds
                if(round > 0) {
                    res = res | (m_x ^ s_x[j]) | (m_z ^ s_z[j]);
                }
                // Update to latest measurement results
                s_x[j] = m_x;
                s_z[j] = m_z;
            }

            // Measure logical-Z to ensure we prepare the logical-|0> state
            bit m_lz = rotated_surface_lz_m();
            // Check parity across rounds
            if(round > 0) {
                res = res | (m_lz ^ s_lz);
            }
            s_lz = m_lz;
        }
        // res will now be 0 if all measurements were the same between rounds
    }

    // HACK!! mwpm_full expects to take Z3 context as first argument in Julia
    // For now, let that pass through ...
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


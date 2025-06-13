// Rotated surface code state preparation using a measurement-free approach in
// arXiv:2303.17211
OPENQASM 3.0;
include "stdgates.inc";

////////////////////////////////////////////////////////////////////////////////
// Conventions
//
// For the rotated surface code, the ucc-ft library assumes the qubits
// are laid out on a grid, where qubits are numbered in row major order,
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
// For arXiv:2303.17211, the layout is mirored, and the qubits are snake-ordered
// to map to a linear physical layout:
// ```
//                 Z
//     q0 •───• q5 ──• q6
//   X    │ Z │  X   |
//     q1 •───• q4 ──• q7
//        │ X │  Z   |    X
//     q2 •───• q3 ──• q8
//          Z
// ```
// The mapping of indicies in this layout to indicies in the ucc-ft layout is:
// [2,5,8,7,4,1,0,3,6] for d=3
// [6,5,0,7,4,1,,8,3,2]  to go from ucc-ft to arXiv:2303.17211 layout.


qubit[9] data;

// Subroutine to do the state preperation circuit in Fig1c
def prepare_state() {

    // Step 0 in Fig1c
    for int i in [0:8] {
        reset data[i];
    }

    for int i in [1:2:8] {
        h data[i];
    }

    // Step 1 in Fig1c
    cx data[1], data[0];
    cx data[3], data[2];
    cx data[5], data[4];
    cx data[7], data[8];
    cx data[3], data[4];
    cx data[5], data[6];

    // Step 2 in Fig1c
    cx data[2], data[1];
    cx data[6], data[7];

    // Looks like fails on cx data[6], data[7] with
    // X and Z error on data[6] and Z error on data[7]
    // [9, 1, 3] code, so t = 1 correctable error

}

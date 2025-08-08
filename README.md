# ucc-ft: Formal Verification for Fault-Tolerant Quantum Circuits

**ucc-ft** is an open-source Python tool for automatically verifying the fault-tolerance of quantum error correction (QEC) circuit implementations. This alpha version wraps and extends the formal verification framework from Chen et al. ([paper](https://arxiv.org/abs/2501.14380v2), [original code](https://zenodo.org/records/15214439)) to work with Python and QASM3.

## Features

- **Automated fault-tolerance verification** for quantum error correction circuit gadgets
- **Standard interfaces** supporting Python and OpenQASM 3 circuit descriptions
- **Symbolic execution** engine for comprehensive error analysis
- **SMT solver integration** for rigorous mathematical verification

## Installation

As alpha software, this is not yet published as a python package. Instead to use, please install [uv](https://docs.astral.sh/uv/getting-started/installation/), which manages dependencies and ensures a reproducible development environment.

```bash
git clone https://github.com/unitaryfoundation/ucc-ft.git
cd ucc-ft
uv sync
```

For development with uv, either prefix commands with `uv run` or activate the virtual environment:
```bash
source .venv/bin/activate
```

For more details on using uv, refer to its [documentation](https://docs.astral.sh/uv/) or [this tutorial](https://realpython.com/python-uv/).

## Quick Start

Here's how to verify the (trivially true) fault tolerance of a transversal CNOT gate for the rotated surface code:

```python
from ucc_ft.codes import RotatedSurfaceCode
from ucc_ft.checker import ft_check

# Define a distance-3 rotated surface code
d3_surface_code = RotatedSurfaceCode(3)

# OpenQASM 3 circuit for transversal CNOT
cnot_circuit = """
OPENQASM 3.0;
include "stdgates.inc";
const uint d = 3;
const uint data_size = d * d;
qubit[data_size] state1;
qubit[data_size] state2;

def logical_CNOT() {
    for int i in [0:(data_size-1)] {
        cx state1[i], state2[i];
    }
}
"""

# Verify fault-tolerance
ft_res = ft_check(d3_surface_code, cnot_circuit, "logical_CNOT", "gate")
print(ft_res)
```

Check out more interesting [./examples here](./examples/).

## How It Works

ucc-ft uses **quantum symbolic execution** to analyze all possible error propagation paths in your quantum circuits:

1. **Circuit Analysis**: Parses OpenQASM 3 circuits and quantum error correction code definitions
2. **Symbolic Execution**: Tracks error propagation through the circuit using symbolic stabilizer states
3. **Verification**: Uses SMT solvers to mathematically prove whether fault-tolerance conditions are satisfied
4. **Results**: Returns verification status and error details when faults are detected

## Dependencies

- **Python dependencies**: Managed automatically via `uv sync`
- **Bitwuzla SMT solver**: Must be installed separately and available on your PATH. Follow the [installation instructions](https://bitwuzla.github.io/docs/index.html) to build from source.

## Acknowledgement

This tool wraps the initial Julia implementation of the core symbolic execution and verification engine in

```bibtex
@article{chen2025verifying,
  title={Verifying Fault-Tolerance of Quantum Error Correction Codes},
  author={Chen, Kean and Liu, Yuxiang and Li, Gushu},
  journal={arXiv preprint arXiv:2501.14380v2},
  year={2025}
}
```

Over time, will be looking to add additional features and capabilities.

## Status

This is a prototype/alpha version under active development. We welcome feedback, bug reports, and contributions!


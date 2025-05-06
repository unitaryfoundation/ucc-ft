# ucc-ft
This repository is a work-in-progress protoype of "Verifying Fault-Tolerance of Quantum Error Correction Codes" by Chen et al.[^kean] ([paper](https://arxiv.org/abs/2501.14380v2), [code](https://zenodo.org/records/15214439)). It follows the approach outlined in [this](https://github.com/unitaryfoundation/ucc/discussions/344) Github discussion, starting from the Julia code linked above, and incrementally translating to python code.

## Running/developing
Please install [uv](https://docs.astral.sh/uv/getting-started/installation/), which is used to managed dependencies and ensure a reproducible development enviroment. Once uv is installed, setup your development environment via

```bash
git clone https://github.com/unitaryfoundation/ucc.git
cd ucc
uv sync
```

For development with uv, we assume you either prefix each command with ``uv run``, or
you first activate the [uv managed virtual environment](https://docs.astral.sh/uv/pip/environments/#using-a-virtual-environment) by running ``source .venv/bin/activate`` in your shell.

For more details on using uv, refer to its [documentation](https://docs.astral.sh/uv/) or [this tutorial](https://realpython.com/python-uv/).

### Bitwuzla Dependencies
In addition to the python dependencies, the code assumes [Bitwuzla](https://bitwuzla.github.io/docs/index.html) was installed and is available on your path. If not, follow the instructions at that link to biuld from source.


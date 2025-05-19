from .julia import jl
from stim import PauliString, Tableau
import numpy as np
import math
import openqasm3
import openqasm3.ast as ast
from openqasm3.printer import Printer, PrinterState
from typing import List, Sequence, Any
import io
from dataclasses import dataclass


class QProgVisitor(Printer):
    """
    Visitor that converts the QASM3 program to an @qprog statement
    that works with QuantumSE.jl in Julia. This is very basic and makes
    specific assumptions about the QASM3 program, raising an error if
    there are any issues.

    For now, we assume the QASM3 program defines a single subroutine, that takes
    two qubit registers as arugments:

      def func(qubit[size] state, qubit ancilla) { BODY }

    The first argument is the state register and corresponds to the physical qubits
    that encode a logical qubit. The ancilla is additional qubits that are used to
    implement the application.

    TODO: Support multiple logical input qubits; larger sized ancilla etc.

    """

    def __init__(self, stream: io.TextIOBase):
        super().__init__(stream)
        self.func_name = None
        self.global_decls = set()

    def visit_Include(self, node: ast.Include, context: PrinterState) -> None:
        # Ignore the include statement when converting to Julia
        pass

    def visit_Program(self, node: ast.Program, context: PrinterState) -> None:
        for statement in node.statements:
            self.visit(statement, context)

    def _visit_statement_list(
        self,
        nodes: List[ast.Statement],
        context: PrinterState,
        prefix: str = "",
    ) -> None:
        # Same as as the base class, but doesn't print '{' and '}' which are not
        # supported in Julia.
        self.stream.write(prefix)
        self._end_line(context)
        with context.increase_scope():
            for statement in nodes:
                self.visit(statement, context)
        self._start_line(context)

    def visit_IntType(self, node: ast.IntType, context: PrinterState) -> None:
        # Julia capitalizes the first letter of the type
        self.stream.write("Int")
        if node.size is not None:
            # This won't work for arbitrary sizes; but can handle later if needed
            self.visit(node.size, context)

    def visit_UintType(self, node: ast.UintType, context: PrinterState) -> None:
        # Julia capitalizes the first letter of the type
        self.stream.write("UInt")
        if node.size is not None:
            # This won't work for arbitrary sizes; but can handle later if needed
            self.visit(node.size, context)

    def visit_BoolType(self, node: ast.BoolType, context: PrinterState) -> None:
        # Julia capitalizes the first letter of the type
        self.stream.write("Bool")

    def visit_BitType(self, node: ast.BitType, context: PrinterState) -> None:
        # treat bits directly as Bool's
        self.stream.write("Bool")
        if node.size is not None:
            raise ValueError("Bit type does not support size")

    def _visit_bit_init(self, node: ast.Expression, context: PrinterState) -> None:
        # Convert bit initialization expressions to Z3 bitvector expressions for us
        # in symbolic execution of @qprog
        if isinstance(node, ast.IntegerLiteral):
            self.stream.write(f"bv_val(ctx, {node.value}, 1)")
        elif isinstance(node, ast.QuantumMeasurement):
            self.visit(node)
        else:
            raise ValueError(
                f"Unsupported bit initialization expression type: {type(node)}"
            )

    def visit_ClassicalDeclaration(
        self, node: ast.ClassicalDeclaration, context: PrinterState
    ) -> None:
        self._start_line(context)
        self.visit(node.identifier, context)
        # No type for Julia
        # self.stream.write("::")
        # self.visit(node.type)
        if node.init_expression is not None:
            self.stream.write(" = ")
            # For bit types, inject the Z3 bittype
            if isinstance(node.type, ast.BitType):
                self._visit_bit_init(node.init_expression, context)
            else:
                self.visit(node.init_expression, context)

        self._end_statement(context)

    def visit_ConstantDeclaration(
        self, node: ast.ConstantDeclaration, context: PrinterState
    ) -> None:
        self._start_line(context)
        self.stream.write("const ")
        self.visit(node.identifier, context)
        # No type for Julia
        # self.stream.write("::")
        # self.visit(node.type)
        self.stream.write(" = ")
        # For bit types, inject the Z3 bittype
        if isinstance(node.type, ast.BitType):
            self._visit_bit_init(node.init_expression, context)
        else:
            self.visit(node.init_expression, context)
        self._end_statement(context)

        # store global constant declarations
        if context.current_indent == 0:
            self.global_decls.add(node.identifier.name)

    def _visit_sequence_and_add_one(
        self,
        nodes: Sequence[ast.QASMNode],
        context: PrinterState,
        *,
        start: str = "",
        end: str = "",
        separator: str,
    ) -> None:
        if start:
            self.stream.write(start)
        for node in nodes[:-1]:
            self.stream.write("( ")
            self.visit(node, context)
            self.stream.write(" ) + 1 ")
            self.stream.write(separator)
        if nodes:
            self.stream.write("( ")
            self.visit(nodes[-1], context)
            self.stream.write(" ) + 1 ")
        if end:
            self.stream.write(end)

    def visit_IndexedIdentifier(
        self, node: ast.IndexedIdentifier, context: PrinterState
    ) -> None:
        # Julia uses 1-based indexing, but QASM is 0-based.
        # We need to convert the index to 1-based for Julia.

        self.visit(node.name, context)
        for index in node.indices:
            self.stream.write("[")
            if isinstance(index, ast.DiscreteSet):
                raise ValueError("DiscreteSet indexing not supported in Julia")
            else:
                self._visit_sequence_and_add_one(index, context, separator=", ")
            self.stream.write("]")

    def visit_QuantumReset(self, node: ast.QuantumReset, context: PrinterState) -> None:
        self._start_line(context)
        self.stream.write("INIT(")
        self.visit(node.qubits, context)
        self.stream.write(")")
        self._end_statement(context)

    def visit_QuantumGate(self, node: ast.QuantumGate, context: PrinterState) -> None:
        self._start_line(context)
        self._visit_gate_identifier(node.name, context)
        if node.arguments:
            ## This may not work for all gates; can revist later if needed
            self._visit_sequence(
                node.arguments, context, start="(", end=")", separator=", "
            )
        self._visit_sequence(node.qubits, context, start="(", end=")", separator=", ")
        self._end_statement(context)

    def _visit_gate_identifier(
        self, node: ast.Identifier, context: PrinterState
    ) -> None:
        # Add gates as needed to map to the types in QuantumSE.jl/SymbolicStabilizer.jl
        translation = {"h": "H", "cx": "CNOT"}
        if node.name in translation:
            self.stream.write(translation[node.name])
        else:
            raise ValueError(f"Unknown gate identifier: {node.name}")

    def visit_RangeDefinition(
        self, node: ast.RangeDefinition, context: PrinterState
    ) -> None:
        # Add parens around each range element
        if node.start is not None:
            self.stream.write("(")
            self.visit(node.start, context)
            self.stream.write(")")
        self.stream.write(":")
        if node.step is not None:
            self.stream.write("(")
            self.visit(node.step, context)
            self.stream.write(")")
            self.stream.write(":")
        if node.end is not None:
            self.stream.write("(")
            self.visit(node.end, context)
            self.stream.write(")")

    def visit_ForInLoop(self, node: ast.ForInLoop, context: PrinterState) -> None:
        self._start_line(context)
        self.stream.write("for ")
        # no type for julia
        # self.visit(node.type)
        # self.stream.write(" ")
        self.visit(node.identifier, context)
        self.stream.write(" in ")
        if isinstance(node.set_declaration, ast.RangeDefinition):
            self.visit(node.set_declaration, context)
        else:
            self.visit(node.set_declaration, context)
        self._visit_statement_list(node.block, context, prefix=" ")
        self.stream.write("end")
        self._end_line(context)

    def visit_QuantumMeasurement(
        self, node: ast.QuantumMeasurement, context: PrinterState
    ) -> None:
        self.stream.write("DestructiveM(")
        self.visit(node.qubit, context)
        self.stream.write(")")

    def visit_SubroutineDefinition(
        self, node: ast.SubroutineDefinition, context: PrinterState
    ) -> None:
        # Here is where we are heavily restricting the QASM3 program format
        # This subroutine is meant to be the main circuit we are checking FT on.
        # We assume it has a specific format:
        #  def func(qubit[size] state, qubit ancilla) { BODY }
        # where state is the logical qubit register and ancilla is any ancilla
        # needed for measurements or other operations.
        if len(node.arguments) != 2:
            raise ValueError("Subroutine must have exactly two arguments")
        if not isinstance(node.arguments[0], ast.QuantumArgument):
            raise ValueError("First argument must be a quantum argument")
        if not isinstance(node.arguments[1], ast.QuantumArgument):
            raise ValueError("Second argument must be a quantum argument")
        if node.arguments[0].name.name != "state":
            raise ValueError("First argument must be named 'state'")
        if node.arguments[1].name.name != "ancilla":
            raise ValueError("Second argument must be named 'ancilla'")

        self.func_name = node.name.name
        self._start_line(context)
        self.stream.write("@qprog ")

        self.visit(node.name, context)
        self.stream.write(" () begin")
        self._end_line(context)

        ## Declare the arrays for state and ancilla
        ## The qprog code has arrays of integer indices for the qubits here
        self.stream.write(f"  state = [i for i in 1:{node.arguments[0].size.name}]")
        self._end_line(context)
        self.stream.write("  ancilla = length(state) + 1")
        self._end_line(context)

        self._visit_statement_list(node.body, context)
        self.stream.write("end")
        self._end_line(context)

    def visit_WhileLoop(self, node: ast.WhileLoop, context: PrinterState) -> None:
        # converts a QASM while loop to a @qprog repeat-until loop
        # this requries inverting the while loop condition

        self._start_line(context)
        self.stream.write("@repeat begin")
        self._visit_statement_list(node.block, context, prefix=" ")

        self.stream.write("end :until ")

        # Going from while(a) { BODY } to repeat { BODY } until !a
        # Restrict the while loop condition to be a binary expression
        # with LHS an identifier, RHS a literal and operator == or !=

        # Aside (it would be nice to just put ! in front of the entire expression but the
        # overloaded Z3 binary vectors don't seem to support that)

        if not isinstance(node.while_condition, ast.BinaryExpression):
            raise ValueError("While condition must be a unary expression")
        if not isinstance(node.while_condition.lhs, ast.Identifier):
            raise ValueError("While condition LHS must be a variable")
        if not isinstance(node.while_condition.rhs, ast.IntegerLiteral):
            raise ValueError("While condition RHS must be a literal")

        def to_e(x):
            return getattr(ast.BinaryOperator, x)

        op_map = {to_e("=="): to_e("!="), to_e("!="): to_e("==")}

        if node.while_condition.op not in op_map:
            raise ValueError(
                "While condition must be a binary expression with == or !="
            )

        target_op = op_map[node.while_condition.op]
        literal = node.while_condition.rhs.value
        self.stream.write(
            f"({node.while_condition.lhs.name} {target_op.name} bv_val(ctx,{literal},1))"
        )
        self._end_line(context)


def to_julia_tableau_fmt(
    tableau: Tableau,
) -> np.ndarray:
    """
    Converts a Stim stabilizer Tableau to the specified Julia boolean matrix format
    expected by the QuantumSE.jl package.

    The Julia format represents stabilizers as a boolean matrix where rows are
    generators and columns are [X_part | Z_part].

    The Stim format returned by tableau.to_numpy() follows the Aaronson-Gottesman
    format in https://arxiv.org/pdf/quant-ph/0406196.

    """
    _x2x, _x2z, z2x, z2z, _x_signs, _z_signs = tableau.to_numpy()

    num_qubits = z2x.shape[0]

    julia_stabilizer = np.zeros((num_qubits, 2 * num_qubits), dtype=bool)
    # We only care how the stabilizer generators for Z (versus entire Tableu)
    julia_stabilizer[:, 0:num_qubits] = z2x
    julia_stabilizer[:, num_qubits : 2 * num_qubits] = z2z

    return julia_stabilizer


@dataclass
class QProgAndContext:
    src: str
    qprog: Any
    global_decls: dict[str, Any]


def qasm_to_qprog(circuit: str) -> QProgAndContext:
    """Convert the circuit to a qprog object, evaluating and returning a handle to
    the julia function object.
    """

    jl.seval("using QuantumSE")
    res = openqasm3.parse(circuit)
    buff = io.StringIO()
    visitor = QProgVisitor(buff)
    visitor.visit(res)
    jl.seval(buff.getvalue())

    return QProgAndContext(
        src=buff.getvalue(),
        qprog=getattr(
            jl, visitor.func_name
        )(),  # <-- evaluate the qprog to get a handle
        global_decls={s: getattr(jl, s) for s in visitor.global_decls},
    )


def julia_source_to_qprog(
    src: str, func_name: str, decls: List[str]
) -> QProgAndContext:
    """Convert the circuit to a qprog object, evaluating and returning a handle to
    the julia function object.

    Args:
        src (str): The Julia source code to be evaluated.
        func_name (str): The name of the function handle that is the main QProg entry
        decls (List[str]): A list of declarations to be included in the global context.
    """
    jl.seval(src)
    return QProgAndContext(
        src=src,
        qprog=getattr(jl, func_name),
        global_decls={s: getattr(jl, s) for s in decls},
    )


def _bits_needed(j):
    return int(math.ceil(math.log2(j + 1)))


#### idealized FT check interface; will get there progressively
def _make_symbolic_state(num_qubits, ctx, num_ancilla, stabilizers, phases):
    return jl.from_stabilizer_py(
        num_qubits,
        to_julia_tableau_fmt(Tableau.from_stabilizers(stabilizers)),
        phases,
        ctx,
        num_ancilla,
    )


def _prepare_states(code, ctx, num_ancilla):
    # Prepare state is physical z stabilizers (not code stabilizers!!!)
    stabilizers = code.physical_z_stabilizers()
    phases = [jl.bv_val(ctx, 0, 1) for _ in range(len(stabilizers))]
    return [
        _make_symbolic_state(code.num_qubits, ctx, num_ancilla, stabilizers, phases)
    ]


def _measure_states(code, ctx, num_ancilla):
    # Measurement states are the symbolic logical-Z states
    stabilizers = code.stabilizers()
    phases = [jl.bv_val(ctx, 0, 1) for _ in range(len(stabilizers))]
    stabilizers = stabilizers + [code.logical_z()]
    phases = phases + [jl.bv_const(ctx, "lz", 1)]
    return [
        _make_symbolic_state(code.num_qubits, ctx, num_ancilla, stabilizers, phases)
    ]


def _correct_states(code, ctx, num_ancilla):
    # Error correction states are +1 logical-X and the symbolic logical-Z states
    # |+> state
    stabilizers_x = code.stabilizers() + [code.logical_x()]
    phases_x = [jl.bv_val(ctx, 0, 1) for _ in range(len(stabilizers_x))]
    # |0/1> state
    stabilizers_z = code.stabilizers() + [code.logical_z()]
    phases_z = [jl.bv_val(ctx, 0, 1) for _ in range(len(code.stabilizers()))] + [
        jl.bv_const(ctx, "lz", 1)
    ]
    return [
        _make_symbolic_state(
            code.num_qubits, ctx, num_ancilla, stabilizers_x, phases_x
        ),
        _make_symbolic_state(
            code.num_qubits, ctx, num_ancilla, stabilizers_z, phases_z
        ),
    ]


def _gate2_states(code, ctx, num_ancilla):
    # Error correction states are +1 logical-X and the symbolic logical-Z states
    # ... but for both qubits
    num_qubits = code.num_qubits
    code_stabilizers_1q = code.stabilizers()
    logical_z_1q = code.logical_z()
    logical_x_1q = code.logical_x()
    # PauliString.operator+ is tensor product
    code_stabilizers_2q = [
        PauliString("I" * num_qubits) + s for s in code.stabilizers()
    ]
    logical_z_2q = PauliString("I" * num_qubits) + code.logical_z()
    logical_x_2q = PauliString("I" * num_qubits) + code.logical_x()
    # LX state
    stabilizers_lx = (
        code_stabilizers_1q + [logical_x_1q] + code_stabilizers_2q + [logical_x_2q]
    )
    phases_lx = [jl.bv_val(ctx, 0, 1) for _ in range(len(stabilizers_lx))]
    # LZ state
    stabilizers_lz = (
        code_stabilizers_1q + [logical_z_1q] + code_stabilizers_2q + [logical_z_2q]
    )
    phases_lz = (
        [jl.bv_val(ctx, 0, 1) for _ in range(len(code_stabilizers_1q))]
        + [jl.bv_const(ctx, "lz1", 1)]
        + [jl.bv_val(ctx, 0, 1) for _ in range(len(code_stabilizers_1q))]
        + [jl.bv_const(ctx, "lz2", 1)]
    )
    return [
        _make_symbolic_state(
            2 * code.num_qubits, ctx, num_ancilla, stabilizers_lx, phases_lx
        ),
        _make_symbolic_state(
            2 * code.num_qubits, ctx, num_ancilla, stabilizers_lz, phases_lz
        ),
    ]


def error_free_symbolic_output(
    code, symbolic_input_state, gadget_type: str, ctx, num_ancilla
):
    match gadget_type:
        case "prepare":
            # logical |0> (no symbols since specific state)
            stabilizers = code.stabilizers() + [code.logical_prep_stabilizer()]
            phases = [jl.bv_val(ctx, 0, 1) for _ in range(len(stabilizers))]
            return _make_symbolic_state(
                code.num_qubits, ctx, num_ancilla, stabilizers, phases
            )
        case "measurement" | "decoder":
            # same as input, but make a copy
            return jl.copy(symbolic_input_state)
            pass
        case "gate":
            # result of running the gate on the input (error free)
            ## TODO: HOW??
            # Assume CNOT for now (yucky)
            res = jl.copy(symbolic_input_state)
            d = code.d
            for i in range(1, d * d + 1):
                jl.CNOT(res, i, i + d * d)
            return res
        case _:
            raise ValueError(f"Invalid gadget type: {gadget_type}")
    pass


def ft_check_ideal(
    code,
    func_and_context: QProgAndContext,
    gadget_type: str,
    NERRS: int = 12,  # TODO: Can this be inferred -- basically log2 number of maximum labeled errors (so how many bits to track it all))
):
    """
    Check if the given circuit is fault tolerant for the given code and gadget type.
    """
    jl.seval("using Z3")
    jl.seval("using QuantumSE")
    ctx = jl.Context()

    num_ancilla = code.d * code.d - 1  # TODO: infer from circuit?

    state_builders = {
        "prepare": _prepare_states,
        "measurement": _measure_states,
        "decoder": _correct_states,
        "gate": _gate2_states,
    }
    if gadget_type not in state_builders:
        raise ValueError(f"Invalid gadget type: {gadget_type}")
    symbolic_input_states = state_builders[gadget_type](code, ctx, num_ancilla)

    for symbolic_input_state in symbolic_input_states:
        # Reset counter used inside julia code for error names
        jl.clearerrcnt()

        symbolic_target_state = error_free_symbolic_output(
            code, symbolic_input_state, gadget_type, ctx, num_ancilla
        )

        cstate = jl.make_cstate(
            {"ctx": ctx, "d": code.d} | func_and_context.global_decls
        )
        num_errors = (code.d - 1) // 2

        num_main_qubits = code.num_qubits
        num_blocks = 1
        if gadget_type == "gate":
            num_main_qubits = 2 * code.num_qubits
            num_blocks = 2

        b_num_main_qubits = _bits_needed(num_main_qubits)

        # inject start state errors
        nerrs_input = jl.bv_val(ctx, 0, b_num_main_qubits)
        if gadget_type != "prepare":
            nerrs_input = jl.inject_errors(
                symbolic_input_state,
                num_main_qubits,
                ctx,
                nerrs_input,
                b_num_main_qubits,
            )

        cfg1 = jl.SymConfig(func_and_context.qprog, cstate, symbolic_input_state, NERRS)

        # Generate configurations and check_FT
        cfgs1 = jl.QuantSymEx(cfg1)
        for cfg in cfgs1:
            # Need to grab these for measurement case
            # TODO: Was hard coded to match surface code # meas_result=cfg.σ[:final_res], meas_gt=lz
            meas_result = (
                getattr(cfg, "σ").get(jl.seval(":final_res"))
                if gadget_type == "measurement"
                else None
            )
            # In Z3, same name and type are the same, so don't need to keep a handle around from when state was created
            meas_gt = (
                jl.bv_const(ctx, "lz", 1) if gadget_type == "measurement" else None
            )

            if not jl.check_FT_py(
                cfg,
                symbolic_target_state,
                num_errors,
                nerrs_input,
                gadget_type,
                num_blocks=num_blocks,
                meas_result=meas_result,
                meas_gt=meas_gt,
            ):
                return False
    return True

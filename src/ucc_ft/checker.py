from .julia import jl
import juliacall as jc
from stim import PauliString, Tableau
import numpy as np
import math
import openqasm3
import openqasm3.ast as ast
from openqasm3.printer import Printer, PrinterState
from openqasm3 import properties
from typing import List, Sequence, Any, Optional
import re
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
        self.stream.write("__qubit_count = 0\n")
        self.func_names = set()
        self.global_decls = set()
        self.extern_decls = set()

    def visit_Include(self, node: ast.Include, context: PrinterState) -> None:
        # Ignore the include statement when converting to Julia
        pass

    def visit_Program(self, node: ast.Program, context: PrinterState) -> None:
        for statement in node.statements:
            self.visit(statement, context)

    def _inject_source_line(self, line: int, context: PrinterState) -> None:
        # Inject a comment with the source line number for error reporting
        self._start_line(context)
        self.stream.write(f"set_source_line({line})")
        self._end_statement(context)

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
        elif isinstance(node, ast.FunctionCall):
            self.visit(node)
        else:
            raise ValueError(
                f"Unsupported bit initialization expression type: {type(node)}"
            )

    def visit_ClassicalDeclaration(
        self, node: ast.ClassicalDeclaration, context: PrinterState
    ) -> None:
        self._inject_source_line(node.span.start_line, context)
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
        elif isinstance(node.type, ast.BitType) and node.type.size is not None:
            # If it's an uninitialized array bit type, intialize the Z3 expression holder
            self.stream.write(" = Vector{Z3.Expr}(undef, ")
            self.visit(node.type.size, context)
            self.stream.write(" )")
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

    def visit_IndexExpression(
        self, node: ast.IndexExpression, context: PrinterState
    ) -> None:
        # Julia uses 1-based indexing, but QASM is 0-based.
        # We need to convert the index to 1-based for Julia.

        if properties.precedence(node.collection) < properties.precedence(node):
            self.stream.write("(")
            self.visit(node.collection, context)
            self.stream.write(")")
        else:
            self.visit(node.collection, context)
        self.stream.write("[")
        if isinstance(node.index, ast.DiscreteSet):
            raise ValueError("DiscreteSet indexing not supported in Julia")
        else:
            self._visit_sequence_and_add_one(node.index, context, separator=", ")
        self.stream.write("]")

    def visit_BinaryExpression(
        self, node: ast.BinaryExpression, context: PrinterState
    ) -> None:
        our_precedence = properties.precedence(node)
        self.stream.write("(")
        # All AST nodes that are built into BinaryExpression are currently left associative.
        if properties.precedence(node.lhs) < our_precedence:
            self.stream.write("(")
            self.visit(node.lhs, context)
            self.stream.write(")")
        else:
            self.visit(node.lhs, context)
        ## Translate operator to Julia form
        # ⊻ is julia xor vs ^ for QASM
        # ÷ is julia integer division vs / for QASM (so this might not always work)
        translations = {"^": "⊻", "/": "÷"}
        if node.op.name in translations:
            self.stream.write(f" {translations[node.op.name]} ")
        else:
            self.stream.write(f" {node.op.name} ")

        if properties.precedence(node.rhs) <= our_precedence:
            self.stream.write("(")
            self.visit(node.rhs, context)
            self.stream.write(")")
        else:
            self.visit(node.rhs, context)
        self.stream.write(")")

    def visit_QuantumReset(self, node: ast.QuantumReset, context: PrinterState) -> None:
        self._inject_source_line(node.span.start_line, context)
        self._start_line(context)
        self.stream.write("INIT(")
        self.visit(node.qubits, context)
        self.stream.write(")")
        self._end_statement(context)

    def visit_QuantumGate(self, node: ast.QuantumGate, context: PrinterState) -> None:
        self._inject_source_line(node.span.start_line, context)
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
        translation = {"h": "H", "cx": "CNOT", "z": "Z", "x": "X", "y": "Y", "cz": "CZ"}
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

    def _rewrite_conditional_pauli(
        self, node: ast.BranchingStatement, context: PrinterState
    ) -> None:
        # Convert an if (foo) pauli(qubit) to an sX, sZ, sY statement in julia
        # qprog
        # e.g. if( syndrome_X) X(qubit) to sX(qubit, syndrome_X)
        # This is a more efficient way for FT checker code in julia to handle the
        # measurement

        self._inject_source_line(node.if_block[0].span.start_line, context)
        self._start_line(context)
        gate = node.if_block[0].name.name
        rewrite = {"x": "sX", "y": "sY", "z": "sZ"}
        assert gate in rewrite.keys()
        self.stream.write(f"{rewrite[gate]}(")
        self._visit_sequence(node.if_block[0].qubits, context, separator=", ")
        self.stream.write(" , ")
        self.visit(node.condition, context)
        self.stream.write(")")
        self._end_line(context)

    def visit_BranchingStatement(
        self, node: ast.BranchingStatement, context: PrinterState
    ) -> None:
        if (
            len(node.if_block) == 1
            and isinstance(node.if_block[0], ast.QuantumGate)
            and len(node.else_block) == 0
        ):
            return self._rewrite_conditional_pauli(node, context)

        self._start_line(context)
        self.stream.write("if (")
        self.visit(node.condition, context)
        self.stream.write(")")
        self._visit_statement_list(node.if_block, context, prefix=" ")
        if node.else_block:
            self.stream.write(" else")
            # Special handling to flatten a perfectly nested structure of
            #   if {...} else { if {...} else {...} }
            # into the simpler
            #   if {...} else if {...} else {...}
            # but only if we're allowed to by our options.
            if (
                self.chain_else_if
                and len(node.else_block) == 1
                and isinstance(node.else_block[0], ast.BranchingStatement)
                and not node.annotations
            ):
                context.skip_next_indent = True
                self.visit(node.else_block[0], context)
                # Don't end the line, because the outer-most `if` block will.
            else:
                self._visit_statement_list(node.else_block, context)
                self.stream.write("end")
                self._end_line(context)
        else:
            self.stream.write("end")
            self._end_line(context)

    def visit_QuantumMeasurementStatement(self, node, context):
        self._inject_source_line(node.span.start_line, context)
        return super().visit_QuantumMeasurementStatement(node, context)

    def visit_QuantumMeasurement(
        self, node: ast.QuantumMeasurement, context: PrinterState
    ) -> None:
        self.stream.write("DestructiveM(")
        self.visit(node.qubit, context)
        self.stream.write(")")

    def visit_QubitDeclaration(
        self, node: ast.QubitDeclaration, context: PrinterState
    ) -> None:
        self._start_line(context)
        self.visit(node.qubit, context)
        self.global_decls.add(node.qubit.name)
        self.stream.write(" = ")

        # Qubit's in qprog are just indices, so convert these to
        # arrays, keeping track of the number of qubits defined so far
        if node.size is not None:
            if not isinstance(node.size, ast.IntegerLiteral) and not isinstance(
                node.size, ast.Identifier
            ):
                raise ValueError(
                    "Qubit size must be an integer literal or identifier for translation to @qprog"
                )
            self.stream.write("[i + __qubit_count for i in 1:(")
            self.visit(node.size, context)
            self.stream.write(")]")
            self._end_line(context)
            self.stream.write("__qubit_count += ")
            self.visit(node.size, context)
        else:
            self.stream.write("__qubit_count + 1 ")
            self._end_line(context)
            self.stream.write("__qubit_count += 1")
        # self._end_line(context)

        self._end_statement(context)

    def visit_SubroutineDefinition(
        self, node: ast.SubroutineDefinition, context: PrinterState
    ) -> None:
        # For now, we are assuming all arguments are classical, with any
        # quantum registered defined/used as globals.
        # For classical arguments, we drop the type when converting to Julia @qprog

        self.func_names.add(node.name.name)
        self._start_line(context)
        self.stream.write("@qprog ")

        self.visit(node.name, context)
        self.stream.write(" (")
        for arg in node.arguments:
            if isinstance(arg, ast.QuantumArgument):
                raise ValueError(
                    "Quantum arguments are not supported in Julia @qprog, they must be globals"
                )
            self.stream.write(arg.name.name)
        self.stream.write(" ) begin")
        self._end_line(context)

        self._visit_statement_list(node.body, context)
        self.stream.write("end")
        self._end_line(context)

    def visit_ExternDeclaration(
        self, node: ast.ExternDeclaration, context: PrinterState
    ) -> None:
        # Ignore extern declarations
        self.extern_decls.add(node.name.name)

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


def to_stim_paulis(
    julia_tableau: np.ndarray,
) -> List[PauliString]:
    """
    Converts a Julia boolean matrix format to a Stim stabilizer Tableau.

    The Julia format represents stabilizers as a boolean matrix where rows are
    generators and columns are [X_part | Z_part].

    """
    num_qubits = julia_tableau.shape[1] // 2
    xpart = julia_tableau[:, 0:num_qubits]
    zpart = julia_tableau[:, num_qubits:]

    return [
        PauliString.from_numpy(xs=xpart[i], zs=zpart[i])
        for i in range(julia_tableau.shape[0])
    ]


@dataclass
class QProgContext:
    qprog_src: str
    global_decls: dict[str, Any]

    def get_qprog(self, func_name: str, *args) -> jc.AnyValue:
        """
        Get the handle to the qprog object
        """
        return getattr(jl, func_name)(*args)  # <-- evaluate the qprog to get a handle


def qasm_to_qprog_source(qasm_source: str) -> str:
    """
    Translate the given qasm_source to the equivalent Julia qprog source

    """
    res = openqasm3.parse(qasm_source)
    buff = io.StringIO()
    visitor = QProgVisitor(buff)
    visitor.visit(res)
    return buff.getvalue()


def qasm_to_qprog(qasm_source: str) -> QProgContext:
    """
    Translate the given qasm_source and evaluate it within Julia to
    generate a qprog with context.

    """
    jl.seval("using Z3")
    jl.seval("using QuantumSE")
    res = openqasm3.parse(qasm_source)
    buff = io.StringIO()
    visitor = QProgVisitor(buff)
    visitor.visit(res)
    jl.seval(buff.getvalue())

    return QProgContext(
        qprog_src=buff.getvalue(),
        global_decls={s: getattr(jl, s) for s in visitor.global_decls},
    )


def julia_source_to_qprog(src: str, decls: List[str]) -> QProgContext:
    """Convert the circuit to a qprog object

    Args:
        src (str): The Julia source code to be evaluated.
        decls (List[str]): A list of declarations to be included in the global context.
    """
    jl.seval(src)
    return QProgContext(
        qprog_src=src,
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


@dataclass(frozen=True)
class FTErrorLocation:
    """
    Represents a specific fault-tolerant error event.
    """

    error_type: str
    qubit_index: int
    source_line: Optional[int]

    def __str__(self):
        """
        Creates a user-friendly string representation for printing.
        """
        # Format the location part of the string
        if self.source_line is not None:
            location_str = f"[Line {self.source_line: >3}]"  # Right-align line number for clean output
        else:
            location_str = "[Init]"

        return f"{location_str: <6} {self.error_type}-Pauli error on Qubit {self.qubit_index}"


class FTCheckResult:
    def __init__(self, is_ft: bool, error_cause: List[FTErrorLocation] = None):
        self.is_ft = is_ft
        self.error_cause = error_cause

    def __bool__(self):
        return self.is_ft

    def __str__(self):
        if self.is_ft:
            return "Circuit is fault-tolerant"

        return (
            "Circuit is not fault-tolerant. Fault-tolerant error locations:\n"
            + "\n".join(str(e) for e in self.error_cause)
        )


def parse_smt_errors(output: str) -> list[FTErrorLocation]:
    """
    Parses the SMT model output to find all asserted errors.
    """
    error_locations = []

    # Regex to find all `define-fun` lines that result in `#b1`
    error_name_pattern = re.compile(r"\(define-fun\s+([^\s]+)\s+.*\s+#b1\)")

    # Regex to parse the components of an error name.
    # It handles both cases: with and without the optional `_L<number>` part.
    # - Group 1: ([XZ]) -> The error type 'X' or 'Z'
    # - Group 2: (\d+) -> The qubit index
    # - Group 3: (?:_L(\d+))? -> An optional non-capturing group for the line number.
    #   The inner (\d+) is the part we actually capture. It will be None if not present.
    name_parser_pattern = re.compile(r"symb_([XZ])error_Q(\d+)(?:_L(\d+))?_.*")

    error_names = error_name_pattern.findall(output)

    for name in error_names:
        match = name_parser_pattern.match(name)
        if match:
            error_type, qubit_str, line_str = match.groups()

            error_obj = FTErrorLocation(
                error_type=error_type,
                # Convert to 0-based index to match QASM
                qubit_index=int(qubit_str) - 1,
                source_line=int(line_str) if line_str is not None else None,
            )
            error_locations.append(error_obj)

    return error_locations


def ft_check(
    code, qasm: str, qasm_func: str, gadget_type: str, num_ancilla: int = None
) -> FTCheckResult:
    qprog_context = qasm_to_qprog(qasm)

    return ft_check_ideal(
        code,
        qprog_context.get_qprog(qasm_func),
        qprog_context,
        gadget_type,
        NERRS=12,
        num_ancilla=num_ancilla,
    )


def ft_check_ideal(
    code,
    qprog_handle: jc.AnyValue,
    qprog_context: QProgContext,
    gadget_type: str,
    NERRS: int = 12,  # TODO: Can this be inferred -- basically log2 number of maximum labeled errors (so how many bits to track it all))
    num_ancilla: int = None,  # TODO: Can this be inferred from the circuit (e.g. number of ancilla qubits used in the circuit?
) -> FTCheckResult:
    """
    Check if the given circuit is fault tolerant for the given code and gadget type.
    """
    jl.seval("using Z3")
    jl.seval("using QuantumSE")
    ctx = jl.Context()

    if num_ancilla is None:
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

        cstate = jl.make_cstate({"ctx": ctx, "d": code.d} | qprog_context.global_decls)
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

        cfg1 = jl.SymConfig(qprog_handle, cstate, symbolic_input_state, NERRS)

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

            res, res_str = jl.check_FT_py(
                cfg,
                symbolic_target_state,
                num_errors,
                nerrs_input,
                gadget_type,
                num_blocks=num_blocks,
                meas_result=meas_result,
                meas_gt=meas_gt,
            )
            if not res:
                # Extract the error information and relate back to the line in the code
                err_locations = parse_smt_errors(res_str)
                return FTCheckResult(False, err_locations)
    return FTCheckResult(True)

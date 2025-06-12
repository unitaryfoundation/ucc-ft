using Z3

# Layer for converting between python and Julia
# mostly deals with type conversions


function from_stabilizer_py(num_main_qubits::Integer, Stabilizer::AbstractArray, phases::AbstractVector, ctx::Z3.Context, num_ancilla::Integer=0)
    shape = size(Stabilizer)
    println("shape: ", shape)
    @assert num_main_qubits == shape[1] "Stabilizer must have the same number of rows as num_main_qubits"

    phases = convert(Vector{Z3.Expr}, phases)
    @assert length(phases) == num_main_qubits "Phases must have the same number of elements as num_main_qubits"

    from_stabilizer(num_main_qubits, convert(Matrix{Bool},Stabilizer), phases, ctx, num_ancilla)
end

function make_cstate(d::AbstractDict{Any, Any})
    CState(Symbol(k) => v for (k, v) in d)
end

function check_FT_py(cfg::SymConfig, ρ0::SymStabilizerState, num_errors::Int64, nerrs_input::Z3.Expr, FT_type::String; meas_result=0, meas_gt=0, num_blocks=1)
    check_FT(cfg.ρ, ρ0, cfg.ϕ, cfg.nerrs, cfg.NERRS, num_errors, nerrs_input, FT_type, `bitwuzla -rwl 1`, use_z3=false,
        meas_result = meas_result, meas_gt = meas_gt, num_blocks= num_blocks)
end

function inject_errors(state::SymStabilizerState, num_main_qubits::Int64, ctx::Z3.Context, nerrs_input::Z3.Expr, b_num_main_qubits::Int64)
    for j in 1:num_main_qubits
        nerrs_input += zeroext(ctx, inject_symbolic_error(state, j), b_num_main_qubits)
    end
    return nerrs_input
end

function to_stabilizer_py(state::SymStabilizerState)
    return Matrix{Bool}(output_stabilizer_tableau(state))
end

using QuantumSE
using Z3
using LinearAlgebra
using SimpleGF2

function rotate(d, idx)
    i = (idx-1)÷d+1
    j = (idx-1)%d+1
    return (d-j)*d+i
end

function _zadj(d, idx)
    @assert idx<=(d*d-1)÷2
    if idx <= (d-1)÷2
        return [idx*2-1, idx*2]
    elseif idx > d*(d-1)÷2
        return [2*idx, 2*idx+1]
    else
        i = (idx-1) ÷ ((d-1)÷2)
        j = (((idx-1) % ((d-1)÷2)) * 2) + 1 + (i%2)
        return [(i-1)*d+j, (i-1)*d+j+1, i*d+j, i*d+j+1]
    end
end

function _xadj(d, idx)
    res = _zadj(d, idx)
    return [rotate(d,res[i]) for i in 1:length(res)]
end

function mwpm_full_x(ctx, d::Integer, s)
    mwpm_full(ctx, d, s, true, [])
end

function mwpm_full_z(ctx, d::Integer, s, ml)
    mwpm_full(ctx, d, s, false, [ml])
end

function mwpm_full(ctx, d::Integer, s, isX::Bool, ml)

    num_qubits = d*d
    num_stabilizers = length(s)
    num_logical = length(ml)
    ml_adj = num_logical == 1 ? [[(i-1)*d+(d+1)÷2 for i in 1:d]] : []

    adj = isX ? _xadj : _zadj

    mat=zeros(Bool, num_stabilizers+num_logical, num_qubits)
    for j in 1:num_stabilizers
        for k in adj(d, j)
            mat[j, k] = 1
        end
    end
    for j in 1:num_logical
        for k in ml_adj[j]
            mat[num_stabilizers+j, k] = 1
        end
    end

    mat = _extend_rows(mat)
    mat_GF2 = GF2.(mat)
    inv_mat = Matrix{UInt64}(inv(mat_GF2))

    s = vcat(s,ml)

    part_solution = Vector{Z3.Expr}(undef, num_qubits)
    for j in 1 : num_qubits
        part_solution[j] = simplify(reduce(⊻, [inv_mat[j, jj] & s[jj] for jj in 1:num_stabilizers+num_logical]))
    end

    part_solution
end

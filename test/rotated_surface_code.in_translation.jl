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

function _zadj_for_qasm(d, idx)
    @assert idx<=(d*d-1)÷2
    if idx <= (d-1)÷2
        return [idx*2-1, idx*2, -1, -1]
    elseif idx > d*(d-1)÷2
        return [2*idx, 2*idx+1, -1, -1]
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

function mwpm(ctx, d::Integer, s, s_type)

    # assertion
    phi1 = bool_val(ctx, true)
    adj = s_type == "X" ? _xadj : _zadj
    r_bv = alloc_symb(ctx, "assert_recovery_$(s_type)", len = d*d)
    r = [extract(r_bv, j-1, j-1) for j in 1:d*d]
    for j in 1:(d*d-1)÷2
        phi1 = phi1 & (s[j] == reduce(⊻, r[[adj(d, j)...]]))
    end
    phi2 = (sum( (x -> zeroext(ctx, x, _range2b(d*d))).(r) ) <= bv_val(ctx, (d-1)÷2, _range2b(d*d)))
    phi = not(forall(r_bv, not(phi1 & phi2)))

    # claim
    ϕ₁ = bool_val(ctx, true)
    adj = s_type == "X" ? _xadj : _zadj
    r = [alloc_symb(ctx, "recovery_$(s_type)_$(j)") for j in 1:d*d]
    for j in 1:(d*d-1)÷2
        ϕ₁ = ϕ₁ & (s[j] == reduce(⊻, r[[adj(d, j)...]]))
    end
    ϕ₂ = (sum( (x -> zeroext(ctx, x, _range2b(d*d))).(r) ) <= bv_val(ctx, (d-1)÷2, _range2b(d*d)))

    (r, (simplify(not(phi)) | (ϕ₁ & ϕ₂)), bool_val(ctx, true), bool_val(ctx, true))
end

function mwpm2(ctx, d::Integer, s_x, s_z)

    num_qubits = d*d
    num_logical_qubits = 1

    mat=zeros(Bool, num_qubits-num_logical_qubits, 2*num_qubits)
    for j in 1:(num_qubits-num_logical_qubits)÷2
        for k in _xadj(d, j)
            mat[j, k] = 1
        end
        for k in _zadj(d, j)
            mat[j+(num_qubits-num_logical_qubits)÷2, k+num_qubits] = 1
        end
    end

    mat = vcat(_canonicalize_mat(mat)[1:num_qubits+num_logical_qubits,:], mat)
    mat_GF2 = GF2.(mat)
    inv_mat = Matrix{UInt64}(inv(mat_GF2))
    null_space_mat = Matrix{UInt64}(nullspace(mat_GF2[num_qubits+num_logical_qubits+1:2*num_qubits,:]))

    s = vcat(s_x,s_z)

    part_solution = Vector{Z3.Expr}(undef, 2 * num_qubits)
    for j in 1 : 2 * num_qubits
        part_solution[j] = simplify(reduce(⊻, [inv_mat[j, num_qubits+num_logical_qubits+jj] & s[jj] for jj in 1:num_qubits-num_logical_qubits]))
    end

    #assertion
    coefs = alloc_symb(ctx, "assert_coefs", len = num_qubits+num_logical_qubits)
    coefs_vec = Vector{Z3.Expr}(undef, num_qubits+num_logical_qubits)
    for j in 1 : num_qubits+num_logical_qubits
        coefs_vec[j] = extract(coefs, j - 1, j - 1)
    end

    Pauli_err_vec = Vector{Z3.Expr}(undef, 2 * num_qubits)
    for k in 1 : 2*num_qubits
        Pauli_err_vec[k] = simplify(part_solution[k] ⊻ reduce(⊻, [null_space_mat[k, j] & coefs_vec[j] for j in 1 : num_qubits+num_logical_qubits]))
    end

    weight_condition = _Pauli_weight(ctx, Pauli_err_vec)[1] <= bv_val(ctx, (d-1)÷2, _range2b(d*d))
    phi = not(forall(coefs, not(weight_condition)))

    #claim
    coefs = alloc_symb(ctx, "claim_coefs", len = num_qubits+num_logical_qubits)
    coefs_vec = Vector{Z3.Expr}(undef, num_qubits+num_logical_qubits)
    for j in 1 : num_qubits+num_logical_qubits
        coefs_vec[j] = extract(coefs, j - 1, j - 1)
    end

    Pauli_err_vec = Vector{Z3.Expr}(undef, 2 * num_qubits)
    for k in 1 : 2*num_qubits
        Pauli_err_vec[k] = simplify(part_solution[k] ⊻ reduce(⊻, [null_space_mat[k, j] & coefs_vec[j] for j in 1 : num_qubits+num_logical_qubits]))
    end

    weight_condition = _Pauli_weight(ctx, Pauli_err_vec)[1] <= bv_val(ctx, (d-1)÷2, _range2b(d*d))

    (Pauli_err_vec, (simplify(not(phi)) | weight_condition), bool_val(ctx, true), bool_val(ctx, true))
end

@qprog prepare_cat_orig (cat_qubits, verify_qubit, d) begin

    @repeat begin

        #INITP(cat_qubits[1])
        # More efficient approach that relies on combined init+CNOT for more succinct description
        #for i in 2:length(cat_qubits)
        #   INIT2CNOT12(cat_qubits[1], cat_qubits[i])
        #end

        INIT(cat_qubits[1])
        H(cat_qubits[1])
        for i in 2:length(cat_qubits)
            INIT(cat_qubits[i])
            CNOT(cat_qubits[1], cat_qubits[i])
        end

        verify = generate_cat_verification(d, length(cat_qubits))
        res = Vector{Z3.Expr}(undef, length(verify)+1)
        res[1] = bv_val(ctx, 0, 1)
        for i in 1:length(verify)
            # Invert comments to use more succinct version from the paper
            INIT(verify_qubit)
            CNOT(cat_qubits[verify[i][1]], verify_qubit)
            #INIT2CNOT12(cat_qubits[verify[i][1]], verify_qubit)
            CNOT(cat_qubits[verify[i][2]], verify_qubit)
            res[i+1] = DestructiveM(verify_qubit)
            #res[i+1] = CNOT12DestructiveM2(cat_qubits[verify[i][2]], verify_qubit)
        end

    end :until (reduce(|, res) == bv_val(ctx, 0, 1))

end

@qprog prepare_cat_internal (cat_qubits, verify_qubit, d) begin

    @repeat begin

        INIT(cat_qubits[1])
        H(cat_qubits[1])
        for i in 2:length(cat_qubits)
            INIT(cat_qubits[i])
            CNOT(cat_qubits[1], cat_qubits[i])
        end

        verify = generate_cat_verification(d, length(cat_qubits))
        res = bv_val(ctx, 0, 1)
        for i in 1:length(verify)
            INIT(verify_qubit)
            CNOT(cat_qubits[verify[i][1]], verify_qubit)
            CNOT(cat_qubits[verify[i][2]], verify_qubit)
            tmp = DestructiveM(verify_qubit)
            res = res | tmp
        end

    end :until (res == bv_val(ctx, 0, 1))

end


function rotated_surface_x_s(d::Integer, idx::Integer)
	s = zeros(Bool, 2*d*d)

	for j in _xadj(d, idx)
		s[j] = true
	end

	s
end

@qprog rotated_surface_x_m (d, idx) begin
    b = _xadj(d, idx)

    len_b = length(b)
    cat_qubits = [d*d+i for i in 1:len_b]
    verify_qubit = d*d+len_b+1
    prepare_cat_internal(cat_qubits, verify_qubit, d)
    # More efficient approach just checks outputs given cat preparation gadget
    # already checked
    #CatPreparationMod(cat_qubits)
    res = Vector{Z3.Expr}(undef, len_b)
    for i in 1:len_b
        CNOT(cat_qubits[i], b[i])
    end

    for i in 1:len_b
        res[i] = MX(cat_qubits[i])
    end

    res = reduce(⊻, res)

    res
end

function rotated_surface_z_s(d::Integer, idx::Integer)
	s = zeros(Bool, 2*d*d)

	for j in (_zadj(d, idx) .+ d*d)
		s[j] = true
	end

	s
end

function rotated_surface_lx(d::Integer)
	s = zeros(Bool, 2*d*d)

	@inbounds @simd for i in 1:d
		s[((d-1)÷2)*d+i] = true
	end

	s
end

function rotated_surface_lz(d::Integer)
	s = zeros(Bool, 2*d*d)

	@inbounds @simd for i in 1:d
		s[d*d+(i-1)*d+(d+1)÷2] = true
	end

	s
end

@qprog _rotated_surface_decoder (d) begin

    t = Int(floor((d-1)/2))
    #t = t - 1

    @repeat begin

        s_x = Matrix{Z3.Expr}(undef, (d*d-1)÷2, t+1)
        s_z = Matrix{Z3.Expr}(undef, (d*d-1)÷2, t+1)


        for i in 1:t+1
            for j in 1:(d*d-1)÷2
                #s_x[j,i] = rotated_surface_x_m(d, j)
                #s_z[j,i] = rotated_surface_z_m(d, j)
                b = _xadj(d, j)
                s_x[j,i] = MultiPauliMeasurement(b, ['X' for kk in 1:length(b)])
                b = _zadj(d, j)
                s_z[j,i] = MultiPauliMeasurement(b, ['Z' for kk in 1:length(b)])
            end
        end

        check_eq_x = Vector{Z3.Expr}(undef, t)
        check_eq_z = Vector{Z3.Expr}(undef, t)

        for i in 1:t
            check_eq_x[i] = reduce(|, [(s_x[j,i] ⊻ s_x[j,i+1]) for j in 1:(d*d-1)÷2])
            check_eq_z[i] = reduce(|, [(s_z[j,i] ⊻ s_z[j,i+1]) for j in 1:(d*d-1)÷2])
        end

        check_eq = (reduce(|, check_eq_x) | reduce(|, check_eq_z))

    end :until (check_eq == bv_val(ctx, 0, 1))

    r_x = mwpm(ctx, d, s_x[:, t+1], "X")
    r_z = mwpm(ctx, d, s_z[:, t+1], "Z")
    #r = mwpm2(ctx, d, s_x, s_z)

    for j in 1:d*d
        #sZ(j, r_x[j])
        #sX(j, r_z[j])
        sPauli(j, r_z[j], r_x[j])
        #sPauli(j, r[j+d*d], r[j])
    end

end

@qprog rotated_surface_decoder (d) begin

    _rotated_surface_decoder(d)

    ancilla = [d*d+1, d*d+2, d*d+3, d*d+4, d*d+5]
    for i in 1:5
        INIT(ancilla[i])
    end

end

@qprog rotated_surface_identity (d) begin

    for i in 1:d*d
        Identity(i)
    end

end

@qprog rotated_surface_CNOT (d) begin

    for i in 1:d*d
        CNOT(i, i+d*d)
    end

end

@qprog rotated_surface_lx_m (d) begin
    b = [((d-1)÷2)*d+i for i in 1:d]

    cat_qubits = [d*d+i for i in 1:d]
    verify_qubit = d*d+d+1
    #prepare_cat(cat_qubits, verify_qubit, d)
    CatPreparationMod(cat_qubits)
    res = Vector{Z3.Expr}(undef, d)
    for i in 1:d
        CNOT(cat_qubits[i], b[i])
    end

    for i in 1:d
        res[i] = MX(cat_qubits[i])
    end

    res = reduce(⊻, res)

    res
end

@qprog rotated_surface_lz_m (d) begin
    b = [(i-1)*d+(d+1)÷2 for i in 1:d]

    cat_qubits = [d*d+i for i in 1:d]
    verify_qubit = d*d+d+1
    # More efficient approach just checks outputs given cat preparation gadget
    prepare_cat_internal(cat_qubits, verify_qubit, d)
    #CatPreparationMod(cat_qubits)
    res = Vector{Z3.Expr}(undef, d)
    for i in 1:d
        CZ(cat_qubits[i], b[i])
    end

    for i in 1:d
        res[i] = MX(cat_qubits[i])
    end

    res = reduce(⊻, res)

    res
end


#function mwpm_full(d::Integer, s, s_type, ml, ml_adj)
#
#    # claim
#    ϕ₁ = bool_val(ctx, true)
#    adj = s_type == "X" ? _xadj : _zadj
#    r = [alloc_symb(ctx, "recovery_$(s_type)_$(j)") for j in 1:d*d]
#    for j in 1:(d*d-1)÷2
#        ϕ₁ = ϕ₁ & (s[j] == reduce(⊻, r[[adj(d, j)...]]))
#    end
#    for j in 1:length(ml)
#        ϕ₁ = ϕ₁ & (ml[j] == reduce(⊻, r[[ml_adj[j]...]]))
#    end
#
#    (r, ϕ₁, bool_val(ctx, true), bool_val(ctx, true))
#end

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

@qprog _rotated_surface_prepare_0 (d) begin

    t = Int(floor((d-1)/2))
    #t = t - 1

    for i in 1:d*d
        INIT(i)
    end

    @repeat begin

        s_x = Matrix{Z3.Expr}(undef, (d*d-1)÷2, t+1)
        s_z = Matrix{Z3.Expr}(undef, (d*d+1)÷2, t+1)


        for i in 1:t+1
            for j in 1:(d*d-1)÷2
                s_x[j,i] = rotated_surface_x_m(d, j)
                s_z[j,i] = rotated_surface_z_m(d, j)
                # More efficient MultiPauliMeasurement relies on having already checked thoes gadgets
                #b = _xadj(d, j)
                #s_x[j,i] = MultiPauliMeasurement(b, ['X' for kk in 1:length(b)])
                #b = _zadj(d, j)
                #s_z[j,i] = MultiPauliMeasurement(b, ['Z' for kk in 1:length(b)])
            end
            s_z[(d*d+1)÷2,i] = rotated_surface_lz_m(d)
            #b = [(ii-1)*d+(d+1)÷2 for ii in 1:d]
            #s_z[(d*d+1)÷2,i] = MultiPauliMeasurement(b, ['Z' for kk in 1:length(b)])
        end

        check_eq_x = Vector{Z3.Expr}(undef, t)
        check_eq_z = Vector{Z3.Expr}(undef, t)

        for i in 1:t
            check_eq_x[i] = reduce(|, [(s_x[j,i] ⊻ s_x[j,i+1]) for j in 1:(d*d-1)÷2])
            check_eq_z[i] = reduce(|, [(s_z[j,i] ⊻ s_z[j,i+1]) for j in 1:(d*d+1)÷2])
        end

        check_eq = (reduce(|, check_eq_x) | reduce(|, check_eq_z))

    end :until (check_eq == bv_val(ctx, 0, 1))

    r_x = mwpm_full(ctx, d, s_x[1:(d*d-1)÷2, t+1], true, [])
    r_z = mwpm_full(ctx, d, s_z[1:(d*d-1)÷2, t+1], false, s_z[(d*d+1)÷2:(d*d+1)÷2, t+1])

    for j in 1:d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
        #sPauli(j, r_z[j], r_x[j])
    end

end

@qprog rotated_surface_prepare_0 (d) begin

    _rotated_surface_prepare_0(d)

    ancilla = [d*d+i for i in 1:max(5,d+1)]
    for i in 1:length(ancilla)
        INIT(ancilla[i])
    end

end

function majority(ctx, s)

    # assert length(s) % 2 == 1
    len_s = length(s)

    half = bv_val(ctx, len_s ÷ 2, _range2b(len_s))

    r = alloc_symb(ctx, "majority")

    phi1 = ((sum([zeroext(ctx, s[i], _range2b(len_s)) for i in 1:len_s]) <= half) == (r == _bv_val(ctx, 0)))

    (r, phi1, bool_val(ctx, true), bool_val(ctx, true))

end

@qprog rotated_surface_measurement (d) begin

    t = Int(floor((d-1)/2))

    m_lz = Vector{Z3.Expr}(undef, 2*t+1)

    for i in 1:2*t+1
        #m_lz[i] = rotated_surface_lz_m(d)
        b = [(ii-1)*d+(d+1)÷2 for ii in 1:d]
        m_lz[i] = MultiPauliMeasurement(b, ['Z' for kk in 1:length(b)])
        _rotated_surface_decoder(d)
    end

    final_res = majority(ctx, m_lz)

end

#######

@qprog rotated_surface_syndrome_measurement (d) begin

    s_x = Vector{Z3.Expr}(undef, (d*d-1)÷2)
    s_z = Vector{Z3.Expr}(undef, (d*d-1)÷2)

    for j in 1:(d*d-1)÷2
        s_x[j] = rotated_surface_x_m(d, j)
        s_z[j] = rotated_surface_z_m(d, j)
        #b = _xadj(d, j)
        #s_x[j] = MultiPauliMeasurement(b, ['X' for kk in 1:length(b)])
        #b = _zadj(d, j)
        #s_z[j] = MultiPauliMeasurement(b, ['Z' for kk in 1:length(b)])
    end

    ancilla = [d*d+1, d*d+2, d*d+3, d*d+4, d*d+5]
    for i in 1:5
        INIT(ancilla[i])
    end

end


"""
    embed(basis1[, basis2], operators::Dict)

`operators` is a dictionary `Dict{Vector{Int}, AbstractOperator}`. The integer vector
specifies in which subsystems the corresponding operator is defined.
"""
function embed(basis_l::CompositeBasis, basis_r::CompositeBasis,
               operators::Dict{<:Vector{<:Integer}, T}) where T<:AbstractOperator
    @assert length(basis_l.bases) == length(basis_r.bases)
    N = length(basis_l.bases)::Int # type assertion to help type inference
    if length(operators) == 0
        return identityoperator(T, basis_l, basis_r)
    end
    indices, operator_list = zip(operators...)
    operator_list = [operator_list...;]
    S = mapreduce(eltype, promote_type, operator_list)
    indices_flat = [indices...;]::Vector{Int} # type assertion to help type inference
    start_indices_flat = [i[1] for i in indices]
    complement_indices_flat = Int[i for i=1:N if i ∉ indices_flat]
    operators_flat = AbstractOperator[]
    if all(([minimum(I):maximum(I);]==I)::Bool for I in indices) # type assertion to help type inference
        for i in 1:N
            if i in complement_indices_flat
                push!(operators_flat, identityoperator(T, S, basis_l.bases[i], basis_r.bases[i]))
            elseif i in start_indices_flat
                push!(operators_flat, operator_list[indexin(i, start_indices_flat)[1]])
            end
        end
        return tensor(operators_flat...)
    else
        complement_operators = [identityoperator(T, S, basis_l.bases[i], basis_r.bases[i]) for i in complement_indices_flat]
        op = tensor([operator_list; complement_operators]...)
        perm = sortperm([indices_flat; complement_indices_flat])
        return permutesystems(op, perm)
    end
end
embed(basis_l::CompositeBasis, basis_r::CompositeBasis, operators::Dict{<:Integer, T}; kwargs...) where {T<:AbstractOperator} = embed(basis_l, basis_r, Dict([i]=>op_i for (i, op_i) in operators); kwargs...)
embed(basis::CompositeBasis, operators::Dict{<:Integer, T}; kwargs...) where {T<:AbstractOperator} = embed(basis, basis, operators; kwargs...)
embed(basis::CompositeBasis, operators::Dict{<:Vector{<:Integer}, T}; kwargs...) where {T<:AbstractOperator} = embed(basis, basis, operators; kwargs...)

# The dictionary implementation works for non-DataOperators
embed(basis_l::CompositeBasis, basis_r::CompositeBasis, indices, op::T) where T<:AbstractOperator = embed(basis_l, basis_r, Dict(indices=>op))

embed(basis_l::CompositeBasis, basis_r::CompositeBasis, index::Integer, op::AbstractOperator) = embed(basis_l, basis_r, index, [op])
embed(basis::CompositeBasis, indices, operators::Vector{T}) where {T<:AbstractOperator} = embed(basis, basis, indices, operators)
embed(basis::CompositeBasis, indices, op::AbstractOperator) = embed(basis, basis, indices, op)

"""
    embed(basis::Basis, indices, ops::AbstractOperator)
    
Special case of the `embed` to handle embedding an operator `ops` into a single
basis.
"""
function embed(basis::Basis, indices, ops::AbstractOperator)
    if indices == 1 || indices == (1,) || indices == [1]
        return ops
    else
        throw(ArgumentError("Invalid indices for single basis."))
    end
end

"""
    embed(basis1[, basis2], indices::Vector, operators::Vector)

Tensor product of operators where missing indices are filled up with identity operators.
"""
function embed(basis_l::CompositeBasis, basis_r::CompositeBasis,
               indices, operators::Vector{T}) where T<:AbstractOperator

    @assert check_embed_indices(indices)

    N = length(basis_l.bases)
    @assert length(basis_r.bases) == N
    @assert length(indices) == length(operators)

    # Embed all single-subspace operators.
    idxop_sb = [x for x in zip(indices, operators) if x[1] isa Integer]
    indices_sb = [x[1] for x in idxop_sb]
    ops_sb = [x[2] for x in idxop_sb]

    for (idxsb, opsb) in zip(indices_sb, ops_sb)
        (opsb.basis_l == basis_l.bases[idxsb]) || throw(IncompatibleBases())
        (opsb.basis_r == basis_r.bases[idxsb]) || throw(IncompatibleBases())
    end

    S = length(operators) > 0 ? mapreduce(eltype, promote_type, operators) : Any
    embed_op = tensor([i ∈ indices_sb ? ops_sb[indexin(i, indices_sb)[1]] : identityoperator(T, S, basis_l.bases[i], basis_r.bases[i]) for i=1:N]...)

    # Embed all joint-subspace operators.
    idxop_comp = [x for x in zip(indices, operators) if x[1] isa Array]
    for (idxs, op) in idxop_comp
        embed_op *= embed(basis_l, basis_r, idxs, op)
    end

    return embed_op
end

permutesystems(a::AbstractOperator, perm) = arithmetic_unary_error("Permutations of subsystems", a)

nsubsystems(s::AbstractKet) = nsubsystems(basis(s))
nsubsystems(s::AbstractOperator) = nsubsystems(basis(s))
nsubsystems(b::CompositeBasis) = length(b.bases)
nsubsystems(b::Basis) = 1
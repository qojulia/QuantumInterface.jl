
"""
    embed(basis1[, basis2], operators::Dict)

`operators` is a dictionary `Dict{Vector{Int}, AbstractOperator}`. The integer vector
specifies in which subsystems the corresponding operator is defined.
"""
function embed(bl::CompositeBasis, br::CompositeBasis,
               operators::Dict{<:Vector{<:Integer}, T}) where T<:AbstractOperator
    (nsubsystems(bl) == nsubsystems(br)) || throw(ArgumentError("Must have nsubsystems(bl) == nsubsystems(br) in embed"))
    N = nsubsystems(bl)::Int # type assertion to help type inference
    if length(operators) == 0
        return identityoperator(T, bl, br)
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
                push!(operators_flat, identityoperator(T, S, bl[i], br[i]))
            elseif i in start_indices_flat
                push!(operators_flat, operator_list[indexin(i, start_indices_flat)[1]])
            end
        end
        return tensor(operators_flat...)
    else
        complement_operators = [identityoperator(T, S, bl[i], br[i]) for i in complement_indices_flat]
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
    embed(basis1[, basis2], indices::Vector, operators::Vector)

Tensor product of operators where missing indices are filled up with identity operators.
"""
function embed(bl::CompositeBasis, br::CompositeBasis,
               indices, operators::Vector{T}) where T<:AbstractOperator

    check_embed_indices(indices) || throw(ArgumentError("Must have unique indices in embed"))
    (nsubsystems(bl) == nsubsystems(br)) || throw(ArgumentError("Must have nsubsystems(bl) == nsubsystems(br) in embed"))
    (length(indices) == length(operators)) || throw(ArgumentError("Must have length(indices) == length(operators) in embed"))

    N = nsubsystems(bl)

    # Embed all single-subspace operators.
    idxop_sb = [x for x in zip(indices, operators) if x[1] isa Integer]
    indices_sb = [x[1] for x in idxop_sb]
    ops_sb = [x[2] for x in idxop_sb]

    for (idxsb, opsb) in zip(indices_sb, ops_sb)
        (basis_l(opsb) == bl[idxsb]) || throw(IncompatibleBases())
        (basis_r(opsb) == br[idxsb]) || throw(IncompatibleBases())
    end

    S = length(operators) > 0 ? mapreduce(eltype, promote_type, operators) : Any
    embed_op = tensor([i ∈ indices_sb ? ops_sb[indexin(i, indices_sb)[1]] : identityoperator(T, S, bl[i], br[i]) for i=1:N]...)

    # Embed all joint-subspace operators.
    idxop_comp = [x for x in zip(indices, operators) if x[1] isa Array]
    for (idxs, op) in idxop_comp
        embed_op *= embed(bl, br, idxs, op)
    end

    return embed_op
end

permutesystems(a::AbstractOperator, perm) = arithmetic_unary_error("Permutations of subsystems", a)

"""
    nsubsystems(a)

Return the number of subsystems of a quantum object in its tensor product
decomposition.

See also [`CompositeBasis`](@ref).
"""
nsubsystems(s::StateVector) = nsubsystems(basis(s))
nsubsystems(s::AbstractOperator) = nsubsystems(basis(s))
nsubsystems(b::CompositeBasis) = length(b.bases)
nsubsystems(b::Basis) = 1
nsubsystems(::Nothing) = 1 # TODO Exists because of QuantumSavory; Consider removing this and reworking the functions that depend on it. E.g., a reason to have it when performing a project_traceout measurement on a state that contains only one subsystem

"""
    nsubspaces(a)

Return the number of subspaces of a quantum object in its direct sum
decomposition.

See also [`SumBasis`](@ref).
"""
nsubspaces(b::SumBasis) = length(b.bases)

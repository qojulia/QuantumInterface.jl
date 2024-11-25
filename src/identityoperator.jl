one(x::Union{<:Basis,<:AbstractOperator}) = identityoperator(x)

"""
    identityoperator(a::Basis[, b::Basis])
    identityoperator(::Type{<:AbstractOperator}, a::Basis[, b::Basis])
    identityoperator(::Type{<:Number}, a::Basis[, b::Basis])
    identityoperator(::Type{<:AbstractOperator}, ::Type{<:Number}, a::Basis[, b::Basis])

Return an identityoperator in the given bases. One can optionally specify the container
type which has to a subtype of [`AbstractOperator`](@ref) as well as the number type
to be used in the identity matrix.
"""
identityoperator(::Type{T}, ::Type{S}, b1::Basis, b2::Basis) where {T<:AbstractOperator,S} = throw(ArgumentError("Identity operator not defined for operator type $T."))
identityoperator(::Type{T}, ::Type{S}, b::OperatorBasis) where {T<:AbstractOperator,S} = identityoperator(T,S,b.left,b.right)
identityoperator(::Type{T}, ::Type{S}, b::Basis) where {T<:AbstractOperator,S} = identityoperator(T,S,b,b)
identityoperator(::Type{T}, b::OperatorBasis) where {T<:AbstractOperator} = identityoperator(T,eltype(T),b)
identityoperator(::Type{T}, bases::Basis...) where {T<:AbstractOperator} = identityoperator(T,eltype(T),bases...)
identityoperator(b::Basis) = identityoperator(ComplexF64,b)
identityoperator(op::T) where {T<:AbstractOperator} = identityoperator(T, fullbasis(op))

# Catch case where eltype cannot be inferred from type; this is a bit hacky
identityoperator(::Type{T}, ::Type{Any}, b1::Basis, b2::Basis) where T<:AbstractOperator = identityoperator(T, ComplexF64, b1, b2)

identityoperator(b1::Basis, b2::Basis) = identityoperator(ComplexF64, b1, b2)

"""Prepare the identity superoperator over a given space."""
function identitysuperoperator end

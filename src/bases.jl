##
# GenericBasis, CompositeBasis
##

"""
    length(b::Basis)

Total dimension of the Hilbert space.
"""
Base.length(b::Basis) = prod(b.shape)

"""
    GenericBasis(N)

A general purpose basis of dimension N.

Should only be used rarely since it defeats the purpose of checking that the
bases of state vectors and operators are correct for algebraic operations.
The preferred way is to specify special bases for different systems.
"""
struct GenericBasis{S} <: Basis
    shape::S
end
GenericBasis(N::Integer) = GenericBasis([N])

Base.:(==)(b1::GenericBasis, b2::GenericBasis) = equal_shape(b1.shape, b2.shape)


"""
    CompositeBasis(b1, b2...)

Basis for composite Hilbert spaces.

Stores the subbases in a vector and creates the shape vector directly
from the shape vectors of these subbases. Instead of creating a CompositeBasis
directly `tensor(b1, b2...)` or `b1 ⊗ b2 ⊗ …` can be used.
"""
struct CompositeBasis{S,B} <: Basis
    shape::S
    bases::B
end
CompositeBasis(bases) = CompositeBasis([length(b) for b ∈ bases], bases)
CompositeBasis(bases::Basis...) = CompositeBasis((bases...,))
CompositeBasis(bases::Vector) = CompositeBasis((bases...,))

Base.:(==)(b1::T, b2::T) where T<:CompositeBasis = equal_shape(b1.shape, b2.shape)

"""
    equal_shape(a, b)

Check if two shape vectors are the same.
"""
function equal_shape(a, b)
    if a === b
        return true
    end
    if length(a) != length(b)
        return false
    end
    for i=1:length(a)
        if a[i]!=b[i]
            return false
        end
    end
    return true
end

"""
    equal_bases(a, b)

Check if two subbases vectors are identical.
"""
function equal_bases(a, b)
    if a===b
        return true
    end
    for i=1:length(a)
        if a[i]!=b[i]
            return false
        end
    end
    return true
end

##
# Common bases
##

"""
    FockBasis(N,offset=0)

Basis for a Fock space where `N` specifies a cutoff, i.e. what the highest
included fock state is. Similarly, the `offset` defines the lowest included
fock state (default is 0). Note that the dimension of this basis is `N+1-offset`.
"""
struct FockBasis{T} <: Basis
    shape::Vector{T}
    N::T
    offset::T
    function FockBasis(N::T,offset::T=0) where T
        if N < 0 || offset < 0 || N <= offset
            throw(DimensionMismatch())
        end
        new{T}([N-offset+1], N, offset)
    end
end

Base.:(==)(b1::FockBasis, b2::FockBasis) = (b1.N==b2.N && b1.offset==b2.offset)


"""
    NLevelBasis(N)

Basis for a system consisting of N states.
"""
struct NLevelBasis{T} <: Basis
    shape::Vector{T}
    N::T
    function NLevelBasis(N::T) where T<:Integer
        if N < 1
            throw(DimensionMismatch())
        end
        new{T}([N], N)
    end
end

Base.:(==)(b1::NLevelBasis, b2::NLevelBasis) = b1.N == b2.N


"""
    PauliBasis(num_qubits::Int)

Basis for an N-qubit space where `num_qubits` specifies the number of qubits.
The dimension of the basis is 2²ᴺ.
"""
struct PauliBasis{S,B} <: Basis
    shape::S
    bases::B
    function PauliBasis(num_qubits::T) where {T<:Integer}
        shape = [2 for _ in 1:num_qubits]
        bases = Tuple(SpinBasis(1//2) for _ in 1:num_qubits)
        return new{typeof(shape),typeof(bases)}(shape, bases)
    end
end

Base.:(==)(pb1::PauliBasis, pb2::PauliBasis) = length(pb1.bases) == length(pb2.bases)


"""
    SpinBasis(n)

Basis for spin-n particles.

The basis can be created for arbitrary spinnumbers by using a rational number,
e.g. `SpinBasis(3//2)`. The Pauli operators are defined for all possible
spin numbers.
"""
struct SpinBasis{S,T} <: Basis
    shape::Vector{T}
    spinnumber::Rational{T}
    function SpinBasis{S}(spinnumber::Rational{T}) where {S,T<:Integer}
        n = numerator(spinnumber)
        d = denominator(spinnumber)
        @assert d==2 || d==1
        @assert n >= 0
        N = numerator(spinnumber*2 + 1)
        new{spinnumber,T}([N], spinnumber)
    end
end
SpinBasis(spinnumber::Rational) = SpinBasis{spinnumber}(spinnumber)
SpinBasis(spinnumber) = SpinBasis(convert(Rational{Int}, spinnumber))

Base.:(==)(b1::SpinBasis, b2::SpinBasis) = b1.spinnumber==b2.spinnumber


"""
    SumBasis(b1, b2...)

Similar to [`CompositeBasis`](@ref) but for the [`directsum`](@ref) (⊕)
"""
struct SumBasis{S,B} <: Basis
    shape::S
    bases::B
end
SumBasis(bases) = SumBasis(Int[length(b) for b in bases], bases)
SumBasis(shape, bases::Vector) = (tmp = (bases...,); SumBasis(shape, tmp))
SumBasis(bases::Vector) = SumBasis((bases...,))
SumBasis(bases::Basis...) = SumBasis((bases...,))

Base.:(==)(b1::T, b2::T) where T<:SumBasis = equal_shape(b1.shape, b2.shape)
Base.:(==)(b1::SumBasis, b2::SumBasis) = false
Base.length(b::SumBasis) = sum(b.shape)

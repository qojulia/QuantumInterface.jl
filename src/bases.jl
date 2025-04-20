"""
Abstract base class for all specialized bases.

The Basis class is meant to specify a basis of the Hilbert space of the
studied system. Besides basis specific information all subclasses must
implement a shape variable which indicates the dimension of the used
Hilbert space. For a spin-1/2 Hilbert space this would be the
vector `[2]`. A system composed of two spins would then have a
shape vector `[2 2]`.

Composite systems can be defined with help of the [`CompositeBasis`](@ref)
class.
"""
abstract type Basis end

"""
    length(b::Basis)

Total dimension of the Hilbert space.
"""
Base.length(b::Basis) = prod(b.shape)

"""
    basis(a)

Return the basis of an object.

If it's ambiguous, e.g. if an operator has a different left and right basis,
an [`IncompatibleBases`](@ref) error is thrown.
"""
function basis end


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

tensor(b::Basis) = b

"""
    tensor(x::Basis, y::Basis, z::Basis...)

Create a [`CompositeBasis`](@ref) from the given bases.

Any given CompositeBasis is expanded so that the resulting CompositeBasis never
contains another CompositeBasis.
"""
tensor(b1::Basis, b2::Basis) = CompositeBasis([length(b1); length(b2)], (b1, b2))
tensor(b1::CompositeBasis, b2::CompositeBasis) = CompositeBasis([b1.shape; b2.shape], (b1.bases..., b2.bases...))
function tensor(b1::CompositeBasis, b2::Basis)
    N = length(b1.bases)
    shape = vcat(b1.shape, length(b2))
    bases = (b1.bases..., b2)
    CompositeBasis(shape, bases)
end
function tensor(b1::Basis, b2::CompositeBasis)
    N = length(b2.bases)
    shape = vcat(length(b1), b2.shape)
    bases = (b1, b2.bases...)
    CompositeBasis(shape, bases)
end
tensor(bases::Basis...) = reduce(tensor, bases)

function Base.:^(b::Basis, N::Integer)
    if N < 1
        throw(ArgumentError("Power of a basis is only defined for positive integers."))
    end
    tensor([b for i=1:N]...)
end

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
Exception that should be raised for an illegal algebraic operation.
"""
mutable struct IncompatibleBases <: Exception end

const BASES_CHECK = Ref(true)

"""
    @samebases

Macro to skip checks for same bases. Useful for `*`, `expect` and similar
functions.
"""
macro samebases(ex)
    return quote
        BASES_CHECK.x = false
        local val = $(esc(ex))
        BASES_CHECK.x = true
        val
    end
end

"""
    samebases(a, b)

Test if two objects have the same bases.
"""
samebases(b1::Basis, b2::Basis) = b1==b2
samebases(b1::Tuple{Basis, Basis}, b2::Tuple{Basis, Basis}) = b1==b2 # for checking superoperators

"""
    check_samebases(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects don't have
the same bases.
"""
function check_samebases(b1, b2)
    if BASES_CHECK[] && !samebases(b1, b2)
        throw(IncompatibleBases())
    end
end


"""
    multiplicable(a, b)

Check if two objects are multiplicable.
"""
multiplicable(b1::Basis, b2::Basis) = b1==b2

function multiplicable(b1::CompositeBasis, b2::CompositeBasis)
    if !equal_shape(b1.shape,b2.shape)
        return false
    end
    for i=1:length(b1.shape)
        if !multiplicable(b1.bases[i], b2.bases[i])
            return false
        end
    end
    return true
end

"""
    check_multiplicable(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are
not multiplicable.
"""
function check_multiplicable(b1, b2)
    if BASES_CHECK[] && !multiplicable(b1, b2)
        throw(IncompatibleBases())
    end
end

"""
    reduced(a, indices)

Reduced basis, state or operator on the specified subsystems.

The `indices` argument, which can be a single integer or a vector of integers,
specifies which subsystems are kept. At least one index must be specified.
"""
function reduced(b::CompositeBasis, indices)
    if length(indices)==0
        throw(ArgumentError("At least one subsystem must be specified in reduced."))
    elseif length(indices)==1
        return b.bases[indices[1]]
    else
        return CompositeBasis(b.shape[indices], b.bases[indices])
    end
end

"""
    ptrace(a, indices)

Partial trace of the given basis, state or operator.

The `indices` argument, which can be a single integer or a vector of integers,
specifies which subsystems are traced out. The number of indices has to be
smaller than the number of subsystems, i.e. it is not allowed to perform a
full trace.
"""
function ptrace(b::CompositeBasis, indices)
    J = [i for i in 1:length(b.bases) if i ∉ indices]
    length(J) > 0 || throw(ArgumentError("Tracing over all indices is not allowed in ptrace."))
    reduced(b, J)
end


"""
    permutesystems(a, perm)

Change the ordering of the subsystems of the given object.

For a permutation vector `[2,1,3]` and a given object with basis `[b1, b2, b3]`
this function results in `[b2, b1, b3]`.
"""
function permutesystems(b::CompositeBasis, perm)
    (nsubsystems(b) == length(perm)) || throw(ArgumentError("Must have nsubsystems(b) == length(perm) in permutesystems"))
    isperm(perm) || throw(ArgumentError("Must pass actual permeutation to permutesystems"))
    CompositeBasis(b.shape[perm], b.bases[perm])
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
            throw(ArgumentError("Fock cutoff and offset must be positive and cutoff must be less than offset"))
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
        N > 0 || throw(ArgumentError("N must be greater than 0"))
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
        d==2 || d==1 || throw(ArgumentError("Can only construct integer or half-integer spin basis"))
        n >= 0 || throw(ArgumentError("Can only construct positive spin basis"))
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

==(b1::T, b2::T) where T<:SumBasis = equal_shape(b1.shape, b2.shape)
==(b1::SumBasis, b2::SumBasis) = false
length(b::SumBasis) = sum(b.shape)

"""
    directsum(b1::Basis, b2::Basis)

Construct the [`SumBasis`](@ref) out of two sub-bases.
"""
directsum(b1::Basis, b2::Basis) = SumBasis(Int[length(b1); length(b2)], Basis[b1, b2])
directsum(b::Basis) = b
directsum(b::Basis...) = reduce(directsum, b)
function directsum(b1::SumBasis, b2::Basis)
    shape = [b1.shape;length(b2)]
    bases = [b1.bases...;b2]
    return SumBasis(shape, (bases...,))
end
function directsum(b1::Basis, b2::SumBasis)
    shape = [length(b1);b2.shape]
    bases = [b1;b2.bases...]
    return SumBasis(shape, (bases...,))
end
function directsum(b1::SumBasis, b2::SumBasis)
    shape = [b1.shape;b2.shape]
    bases = [b1.bases...;b2.bases...]
    return SumBasis(shape, (bases...,))
end

embed(b::SumBasis, indices, ops) = embed(b, b, indices, ops)

##
# show methods
##

function show(stream::IO, x::GenericBasis)
    if length(x.shape) == 1
        write(stream, "Basis(dim=$(x.shape[1]))")
    else
        s = replace(string(x.shape), " " => "")
        write(stream, "Basis(shape=$s)")
    end
end

function show(stream::IO, x::CompositeBasis)
    write(stream, "[")
    for i in 1:length(x.bases)
        show(stream, x.bases[i])
        if i != length(x.bases)
            write(stream, " ⊗ ")
        end
    end
    write(stream, "]")
end

function show(stream::IO, x::SpinBasis)
    d = denominator(x.spinnumber)
    n = numerator(x.spinnumber)
    if d == 1
        write(stream, "Spin($n)")
    else
        write(stream, "Spin($n/$d)")
    end
end

function show(stream::IO, x::FockBasis)
    if iszero(x.offset)
        write(stream, "Fock(cutoff=$(x.N))")
    else
        write(stream, "Fock(cutoff=$(x.N), offset=$(x.offset))")
    end
end

function show(stream::IO, x::NLevelBasis)
    write(stream, "NLevel(N=$(x.N))")
end

function show(stream::IO, x::SumBasis)
    write(stream, "[")
    for i in 1:length(x.bases)
        show(stream, x.bases[i])
        if i != length(x.bases)
            write(stream, " ⊕ ")
        end
    end
    write(stream, "]")
end

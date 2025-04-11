"""
Abstract type for all specialized bases of a Hilbert space.

This type specifies an orthonormal basis for the Hilbert space of the given
system. All subtypes must implement `Base.:(==)` and `Base.length`, where the
latter should return the total dimension of the Hilbert space

Composite systems can be defined with help of [`CompositeBasis`](@ref).

All relevant properties of concrete subtypes of `Basis` defined in
`QuantumInterface` should be accessed using their documented functions and
should not assume anything about the internal representation of instances of
these types (i.e. do not access the fields of the structs directly).
"""
abstract type Basis end

"""
    basis(a)

Return the basis of a quantum object.

If it's ambiguous, e.g. if an operator has a different left and right basis, an
[`IncompatibleBases`](@ref) error is thrown.

See [`StateVector`](@ref) and [`AbstractOperator`](@ref)
"""
function basis end

"""
    basis_l(a)

Return the left basis of an operator.
"""
function basis_l end

"""
    basis_r(a)

Return the right basis of an operator.
"""
function basis_r end


"""
    length(b::Basis)

Total dimension of the Hilbert space.
"""
Base.length(b::Basis) = throw(ArgumentError("Base.length() is not defined for $(typeof(b))"))

"""
    size(b::Basis)

A vector containing the local dimensions of each Hilbert space in its tensor
product decomposition into subsystems.

See also [`nsubsystems`](@ref) and [`CompositeBasis`](@ref).
"""
Base.size(b::Basis) = [length(b)]

"""
    getindex(b::Basis)

Get the i'th factor in the tensor product decomposition of the basis into
subsystems.

See also [`nsubsystems`](@ref) and [`CompositeBasis`](@ref).
"""
Base.getindex(b::Basis, i) = i==1 ? b : throw(BoundsError("attempted to access a nonexistent subsystem basis"))

##
# GenericBasis, CompositeBasis, SumBasis
##

"""
    GenericBasis(N)

A general purpose basis of dimension N.

Should only be used rarely since it defeats the purpose of checking that the
bases of state vectors and operators are correct for algebraic operations.
The preferred way is to specify special bases for different systems.
"""
struct GenericBasis{S} <: Basis
    dim::S
end
Base.:(==)(b1::GenericBasis, b2::GenericBasis) = b1.dim == b2.dim
Base.length(b::GenericBasis) = b.dim

"""
    CompositeBasis(b1, b2...)

Basis for composite Hilbert spaces.

Stores the subbases in a vector and creates the shape vector directly from the
dimensions of these subbases. Instead of creating a CompositeBasis directly,
`tensor(b1, b2...)` or `b1 ⊗ b2 ⊗ …` should be used.
"""
struct CompositeBasis{B<:Basis,S<:Integer} <: Basis
    shape::Vector{S}
    bases::Vector{B}
end
CompositeBasis(bases) = CompositeBasis([length(b) for b in bases], bases)
CompositeBasis(bases::Basis...) = CompositeBasis([bases...])
CompositeBasis(bases::Tuple) = CompositeBasis([bases...])

Base.:(==)(b1::CompositeBasis, b2::CompositeBasis) = all(((i, j),) -> i == j, zip(b1.bases, b2.bases))
Base.length(b::CompositeBasis) = prod(b.shape)
Base.size(b::CompositeBasis) = b.shape
Base.getindex(b::CompositeBasis, i) = b.bases[i]

"""
    tensor(x::Basis, y::Basis, z::Basis...)

Create a [`CompositeBasis`](@ref) from the given bases.

Any given CompositeBasis is expanded so that the resulting CompositeBasis never
contains another CompositeBasis.
"""
tensor(b1::Basis, b2::Basis) = CompositeBasis([length(b1), length(b2)], [b1, b2])
tensor(b1::CompositeBasis, b2::CompositeBasis) = CompositeBasis([b1.shape; b2.shape], [b1.bases; b2.bases])
tensor(b1::CompositeBasis, b2::Basis) = CompositeBasis([b1.shape; length(b2)], [b1.bases; b2])
tensor(b1::Basis, b2::CompositeBasis) = CompositeBasis([length(b1); b2.shape], [b1; b2.bases])
tensor(bases::Basis...) = reduce(tensor, bases)
tensor(basis::Basis) = basis

function Base.:^(b::Basis, N::Integer)
    if N < 1
        throw(ArgumentError("Power of a basis is only defined for positive integers."))
    end
    tensor([b for i=1:N]...)
end

"""
    SumBasis(b1, b2...)

Similar to [`CompositeBasis`](@ref) but for the [`directsum`](@ref) (⊕)
"""
struct SumBasis{S<:Integer,B<:Basis} <: Basis
    shape::Vector{S}
    bases::Vector{B}
end
SumBasis(bases) = SumBasis([length(b) for b in bases], bases)
SumBasis(bases::Basis...) = SumBasis([bases...])
SumBasis(bases::Tuple) = SumBasis([bases...])

Base.:(==)(b1::SumBasis, b2::SumBasis) = all(((i, j),) -> i == j, zip(b1.bases, b2.bases))
Base.length(b::SumBasis) = sum(b.shape)
Base.getindex(b::SumBasis, i) = getindex(b.bases, i)

"""
    directsum(b1::Basis, b2::Basis)

Construct the [`SumBasis`](@ref) out of two sub-bases.
"""
directsum(b1::Basis, b2::Basis) = SumBasis([length(b1), length(b2)], [b1, b2])
directsum(b1::SumBasis, b2::SumBasis) = SumBasis([b1.shape, b2.shape], [b1.bases; b2.bases])
directsum(b1::SumBasis, b2::Basis) = SumBasis([b1.shape; length(b2)], [b1.bases; b2])
directsum(b1::Basis, b2::SumBasis) = SumBasis([length(b1); b2.shape], [b1; b2.bases])
directsum(bases::Basis...) = reduce(directsum, bases)
directsum(basis::Basis) = basis

embed(b::SumBasis, indices, ops) = embed(b, b, indices, ops)

##
# Basis checks
##

"""
Exception that should be raised for an illegal algebraic operation.
"""
mutable struct IncompatibleBases <: Exception end

const BASES_CHECK = Ref(true)

"""
    @compatiblebases

Macro to skip checks for compatible bases. Useful for `*`, `expect` and similar
functions.
"""
macro compatiblebases(ex)
    return quote
        BASES_CHECK[] = false
        local val = $(esc(ex))
        BASES_CHECK[] = true
        val
    end
end

"""
    samebases(b1::Basis, b2::Basis)

Test if two bases are the same. Equivalant to `==`. See
[`check_samebases`](@ref).
"""
samebases(b1::Basis, b2::Basis) = b1==b2

"""
    check_samebases(a, b)

Throw an [`IncompatibleBases`](@ref) error if the bases are not the same. See
[`samebases`](@ref).
"""
function check_samebases(b1, b2)
    if BASES_CHECK[] && !samebases(b1, b2)
        throw(IncompatibleBases())
    end
end

"""
    check_addible(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are not addible as
determined by `addible(a, b)`.  Disabled by use of [`@compatiblebases`](@ref)
anywhere further up in the call stack.
"""
function check_addible(a, b)
    if BASES_CHECK[] && !addible(a, b)
        throw(IncompatibleBases())
    end
end

"""
    check_multiplicable(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are not multiplicable
as determined by `multiplicable(a, b)`.  Disabled by use of
[`@compatiblebases`](@ref) anywhere further up in the call stack.
"""
function check_multiplicable(a, b)
    if BASES_CHECK[] && !multiplicable(a, b)
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

_index_complement(b::CompositeBasis, indices) = complement(length(b.bases), indices)
reduced(a, indices) = ptrace(a, _index_complement(basis(a), indices))

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
included fock state is. Similarly, the `offset` defines the lowest included fock
state (default is 0). Note that the dimension of this basis is `N+1-offset`.
The [`cutoff`](@ref) and [`offset`](@ref) functions can be used to obtain the
respective properties of a given `FockBasis`.
"""
struct FockBasis{T<:Integer} <: Basis
    N::T
    offset::T
    function FockBasis(N::T,offset::T=0) where T
        if N < 0 || offset < 0 || N <= offset
            throw(DimensionMismatch())
        end
        new{T}(N, offset)
    end
end

Base.:(==)(b1::FockBasis, b2::FockBasis) = (b1.N==b2.N && b1.offset==b2.offset)
Base.length(b::FockBasis) = b.N - b.offset + 1

"""
    cutoff(b::FockBasis)

Return the fock cutoff of the given fock basis.

See [`FockBasis`](@ref).
"""
cutoff(b::FockBasis) = b.N

"""
    offset(b::FockBasis)

Return the offset of the given fock basis.

See [`FockBasis`](@ref).
"""
offset(b::FockBasis) = b.offset


"""
    NLevelBasis(N)

Basis for a system consisting of N states.
"""
struct NLevelBasis{T<:Integer} <: Basis
    N::T
    function NLevelBasis(N::T) where T
        if N < 1
            throw(DimensionMismatch())
        end
        new{T}(N)
    end
end

Base.:(==)(b1::NLevelBasis, b2::NLevelBasis) = b1.N == b2.N
Base.length(b::NLevelBasis) = b.N

"""
    SpinBasis(n)

Basis for spin-n particles.

The basis can be created for arbitrary spin numbers by using a rational number,
e.g. `SpinBasis(3//2)`. The Pauli operators are defined for all possible spin
numbers. The [`spinnumber`](@ref) function can be used to get the spin number
for a `SpinBasis`.
"""
struct SpinBasis{T<:Integer} <: Basis
    spinnumber::Rational{T}
    function SpinBasis(spinnumber::Rational{T}) where T
        n = numerator(spinnumber)
        d = denominator(spinnumber)
        d==2 || d==1 || throw(ArgumentError("Can only construct integer or half-integer spin basis"))
        n >= 0 || throw(ArgumentError("Can only construct positive spin basis"))
        N = numerator(spinnumber*2 + 1)
        new{T}(spinnumber)
    end
end
SpinBasis(spinnumber) = SpinBasis(convert(Rational{Int}, spinnumber))

Base.:(==)(b1::SpinBasis, b2::SpinBasis) = b1.spinnumber==b2.spinnumber
Base.length(b::SpinBasis) = numerator(b.spinnumber*2 + 1)

"""
    spinnumber(b::SpinBasis)

Return the spin number of the given spin basis.

See [`SpinBasis`](@ref).
"""
spinnumber(b::SpinBasis) = b.spinnumber


##
# Operator Bases
##

"""
    KetBraBasis(BL,BR)

The "Ket-Bra" operator basis is the standard representation for the left and
right bases of superoperators. This basis is formed by "vec'ing" the
outer-product "Ket-Bra" basis for an operator with a left Bra basis and right
Ket basis which practically means flipping the Bra to a Ket. The operator itself
is then represented as a "Super-Bra" in this basis and corresponds to
column-stacking its matrix.
"""
struct KetBraBasis{BL<:Basis, BR<:Basis} <: Basis
    left::BL
    right::BR
end
KetBraBasis(b::Basis) = KetBraBasis(b,b)
basis_l(b::KetBraBasis) = b.left
basis_r(b::KetBraBasis) = b.right
Base.:(==)(b1::KetBraBasis, b2::KetBraBasis) = (b1.left == b2.left && b1.right == b2.right)
Base.length(b::KetBraBasis) = length(b.left)*length(b.right)

struct ChoiBasis{BL<:Basis, BR<:Basis} <: Basis
    ref::BL
    out::BR
end
basis_l(b::ChoiBasis) = b.ref
basis_r(b::ChoiBasis) = b.out
Base.:(==)(b1::ChoiBasis, b2::ChoiBasis) = (b1.ref == b2.ref && b1.out == b2.out)
Base.length(b::ChoiBasis) = length(b.ref)*length(b.out)

"""
    PauliBasis(N)

The standard Pauli operator basis for an `N` qubit space. This consists of
tensor products of the Pauli matrices I, X, Y, Z, in that order for each qubit.
The dimension of the basis is 2²ᴺ.
"""
struct PauliBasis{T<:Integer} <: Basis
    N::T
end
Base.:(==)(b1::PauliBasis, b2::PauliBasis) = b1.N == b2.N
Base.length(b::PauliBasis) = 4^b.N

"""
    HWPauliBasis(N)

The Hesienberg-Weyl Pauli operator basis consisting of the
N represents the underlying Hilbert
space dimension, not the operator basis dimension. For N>2, this representes the
operator basis formed by the generalized Pauli matrices, also called the clock
and shift matrices. The ordering is the usual one: when the index is written in
base-N and thus has only two digits, the least significant bit gives powers of Z
(the clock matrix), and most significant bit gives powers of X (the shfit matrix).
"""
struct HWPauliBasis{T<:Integer} <: Basis
    shape::Vector{T}
end
HWPauliBasis(N::Integer) = HWPauliBasis([N])
Base.:(==)(b1::HWPauliBasis, b2::HWPauliBasis) = b1.shape == b2.shape
Base.length(b::HWPauliBasis) = prod(b.shape)
Base.getindex(b::HWPauliBasis, i) = HWPauliBasis([b.shape[i]])


##
# show methods
##

function show(stream::IO, x::GenericBasis)
    write(stream, "Basis(dim=$(x.dim))")
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

function show(stream::IO, x::KetBraBasis)
    write(stream, "KetBra(left=$(x.left), right=$(x.right))")
end

function show(stream::IO, x::PauliBasis)
    write(stream, "Pauli(N=$(x.N)")
end

function show(stream::IO, x::HWPauliBasis)
    write(stream, "Pauli($(x.shape)")
end

"""
Abstract type for all specialized bases of a Hilbert space.

This type specifies an orthonormal basis for the Hilbert space of the given
system. All subtypes must implement `Base.:(==)` and `dimension`, where the
latter should return the total dimension of the Hilbert space.

Composite systems can be defined with help of [`CompositeBasis`](@ref).
Custom subtypes can also define composite systems by implementing
`Base.length` and `Base.getindex`.

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

Return the number of subsystems of a quantum object in its tensor product
decomposition.

See also [`CompositeBasis`](@ref).
"""
Base.length(b::Basis) = 1

"""
    getindex(b::Basis)

Get the i'th factor in the tensor product decomposition of the basis into
subsystems.

See also [`CompositeBasis`](@ref).
"""
Base.getindex(b::Basis, i) = i==1 ? b : throw(BoundsError(b,i))
Base.firstindex(b::Basis) = 1
Base.lastindex(b::Basis) = length(b)

Base.iterate(b::Basis, state=1) = state > length(b) ? nothing : (b[state], state+1)

"""
    dimension(b::Basis)

Total dimension of the Hilbert space.
"""
dimension(b::Basis) = throw(ArgumentError("dimesion() is not defined for $(typeof(b))"))

"""
    shape(b::Basis)

A vector containing the local dimensions of each Hilbert space in its tensor
product decomposition into subsystems.

See also [`CompositeBasis`](@ref).
"""
shape(b::Basis) = [dimension(b[i]) for i=1:length(b)]

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
dimension(b::GenericBasis) = b.dim

"""
    CompositeBasis(b1, b2...)

Basis for composite Hilbert spaces.

Stores the subbases in a vector and creates the shape vector directly from the
dimensions of these subbases. Instead of creating a CompositeBasis directly,
`tensor(b1, b2...)` or `b1 ⊗ b2 ⊗ …` should be used.
"""
struct CompositeBasis{B<:Basis} <: Basis
    bases::Vector{B}
    shape::Vector{Int}
    lengths::Vector{Int}
    N::Int
    D::Int
    function CompositeBasis(bases::Vector{B}) where B<:Basis
        # to enable this check the the lazy operators in QuantumOpticsBase need to be changed
        #length(bases) > 1 || throw(ArgumentError("CompositeBasis must only be used for composite systems"))
        shape_ = mapreduce(shape, vcat, bases)
        lengths = cumsum(map(length, bases))
        new{B}(bases, shape_, lengths, lengths[end], prod(shape_))
    end
end
CompositeBasis(bases::Basis...) = CompositeBasis([bases...])
CompositeBasis(bases::Tuple) = CompositeBasis([bases...])

Base.:(==)(b1::CompositeBasis, b2::CompositeBasis) = all(((i, j),) -> i == j, zip(b1.bases, b2.bases))
Base.length(b::CompositeBasis) = b.N
function Base.getindex(b::CompositeBasis, i::Integer)
    (i < 1 || i > b.N) && throw(BoundsError(b,i))
    bases_idx = findfirst(l -> i<=l, b.lengths) 
    inner_idx = i - (bases_idx == 1 ? 0 : b.lengths[bases_idx-1])
    b.bases[bases_idx][inner_idx]
end
Base.getindex(b::CompositeBasis, indices) = [b[i] for i in indices]
shape(b::CompositeBasis) = b.shape
dimension(b::CompositeBasis) = b.D

"""
    tensor(x::Basis, y::Basis, z::Basis...)

Create a [`CompositeBasis`](@ref) from the given bases.

Any given CompositeBasis is expanded so that the resulting CompositeBasis never
contains another CompositeBasis.
"""
tensor(b1::Basis, b2::Basis) = CompositeBasis([b1, b2])
tensor(bases::Basis...) = reduce(tensor, bases)
tensor(basis::Basis) = basis

function tensor(b1::CompositeBasis, b2::CompositeBasis)
    if typeof(b1.bases[end]) == typeof(b2.bases[1])
        t = tensor(b1.bases[end], b2.bases[1])
        if !(t isa CompositeBasis)
            return CompositeBasis([b1.bases[1:end-1]; t;  b2.bases[2:end]])
        end
    end
    return CompositeBasis([b1.bases; b2.bases])
end

function tensor(b1::CompositeBasis, b2::Basis)
    if b1.bases[end] isa typeof(b2)
        t = tensor(b1.bases[end], b2)
        if !(t isa CompositeBasis)
            return CompositeBasis([b1.bases[1:end-1]; t])
        end
    end
    return CompositeBasis([b1.bases; b2])
end

function tensor(b1::Basis, b2::CompositeBasis)
    if b2.bases[1] isa typeof(b1)
        t = tensor(b1, b2.bases[1])
        if !(t isa CompositeBasis)
            return CompositeBasis([t; b2[2:end]])
        end
    end
    return CompositeBasis([b1; b2.bases])
end

Base.:^(b::Basis, N::Integer) = tensor_pow(b, N)

"""
    SumBasis(b1, b2...)

Similar to [`CompositeBasis`](@ref) but for the [`directsum`](@ref) (⊕)
"""
struct SumBasis{S<:Integer,B<:Basis} <: Basis
    shape::Vector{S}
    bases::Vector{B}
end
SumBasis(bases) = SumBasis([dimension(b) for b in bases], bases)
SumBasis(bases::Basis...) = SumBasis([bases...])
SumBasis(bases::Tuple) = SumBasis([bases...])

Base.:(==)(b1::SumBasis, b2::SumBasis) = all(((i, j),) -> i == j, zip(b1.bases, b2.bases))
dimension(b::SumBasis) = sum(b.shape)



"""
    nsubspaces(b)

Return the number of subspaces of a [`SumBasis`](@ref) in its direct sum
decomposition.
"""
nsubspaces(b::SumBasis) = length(b.bases)

"""
    subspace(b, i)

Return the basis for the `i`th subspace of of a [`SumBasis`](@ref).
"""
subspace(b::SumBasis, i) = b.bases[i]

"""
    directsum(b1::Basis, b2::Basis)

Construct the [`SumBasis`](@ref) out of two sub-bases.
"""
directsum(b1::Basis, b2::Basis) = SumBasis([dimension(b1), dimension(b2)], [b1, b2])
directsum(b1::SumBasis, b2::SumBasis) = SumBasis([b1.shape, b2.shape], [b1.bases; b2.bases])
directsum(b1::SumBasis, b2::Basis) = SumBasis([b1.shape; dimension(b2)], [b1.bases; b2])
directsum(b1::Basis, b2::SumBasis) = SumBasis([dimension(b1); b2.shape], [b1; b2.bases])
directsum(bases::Basis...) = reduce(directsum, bases)
directsum(basis::Basis) = basis

# TODO: what to do about embed for SumBasis?
#embed(b::SumBasis, indices, ops) = embed(b, b, indices, ops)

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
function reduced(b::Basis, indices)
    if length(indices)==0
        throw(ArgumentError("At least one subsystem must be specified in reduced."))
    elseif length(indices)==1
        return b[indices[1]]
    else
        return tensor(b[indices]...)
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
function ptrace(b::Basis, indices)
    J = [i for i in 1:length(b) if i ∉ indices]
    length(J) > 0 || throw(ArgumentError("Tracing over all indices is not allowed in ptrace."))
    reduced(b, J)
end

_index_complement(b::Basis, indices) = complement(length(b), indices)
reduced(a, indices) = ptrace(a, _index_complement(basis(a), indices))

"""
    permutesystems(a, perm)

Change the ordering of the subsystems of the given object.

For a permutation vector `[2,1,3]` and a given object with basis `[b1, b2, b3]`
this function results in `[b2, b1, b3]`.
"""
function permutesystems(b::Basis, perm)
    (length(b) == length(perm)) || throw(ArgumentError("Must have length(b) == length(perm) in permutesystems"))
    isperm(perm) || throw(ArgumentError("Must pass actual permeutation to permutesystems"))
    tensor(b.bases[perm]...)
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
            throw(ArgumentError("Fock cutoff and offset must be positive and cutoff must be less than offset"))
        end
        new{T}(N, offset)
    end
end

Base.:(==)(b1::FockBasis, b2::FockBasis) = (b1.N==b2.N && b1.offset==b2.offset)
dimension(b::FockBasis) = b.N - b.offset + 1

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
        N > 0 || throw(ArgumentError("N must be greater than 0"))
        new{T}(N)
    end
end

Base.:(==)(b1::NLevelBasis, b2::NLevelBasis) = b1.N == b2.N
dimension(b::NLevelBasis) = b.N

"""
    SpinBasis(n, N=1)

Basis for spin-`n` particles over `N` systems.

The basis can be created for arbitrary spin numbers by using a rational number,
e.g. `SpinBasis(3//2)`. The Pauli operators are defined for all possible spin
numbers. The [`spinnumber`](@ref) function can be used to get the spin number
for a `SpinBasis`.
"""
struct SpinBasis{T<:Integer} <: Basis
    spinnumber::Rational{T}
    D::T
    N::T
    function SpinBasis(spinnumber::Rational{T}, N=1) where T
        n = numerator(spinnumber)
        d = denominator(spinnumber)
        d==2 || d==1 || throw(ArgumentError("Can only construct integer or half-integer spin basis"))
        n >= 0 || throw(ArgumentError("Can only construct positive spin basis"))
        D = numerator(spinnumber*2 + 1)
        new{T}(spinnumber, D, N)
    end
end
SpinBasis(spinnumber) = SpinBasis(convert(Rational{Int}, spinnumber))

Base.:(==)(b1::SpinBasis, b2::SpinBasis) = b1.D==b2.D && b1.N == b2.N
Base.length(b::SpinBasis) = b.N
Base.getindex(b::SpinBasis, i) = SpinBasis(b.spinnumber, length(i))
shape(b::SpinBasis) = fill(b.D, b.N)
dimension(b::SpinBasis) = b.D^b.N
function tensor(b1::SpinBasis, b2::SpinBasis)
    if b1.spinnumber == b2.spinnumber
        return SpinBasis(b1.spinnumber, b1.N+b2.N)
    else
        return CompositeBasis([b1, b2])
    end
end

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
dimension(b::KetBraBasis) = dimension(b.left)*dimension(b.right)
tensor(b1::KetBraBasis, b2::KetBraBasis) = KetBraBasis(tensor(b1.left,b2.left), tensor(b1.right, b2.right))

"""
    ChoiBasis(ref_basis,out_basis)

The Choi basis is used to represent superoperators in the Choi representation
where the `ref_basis` denotes the ancillary reference system with which an input
state will be jointly measured in order to accomplish teleportation simulation
of the channel with the channel's output appearing in the `out_basis` system.
"""
struct ChoiBasis{BL<:Basis, BR<:Basis} <: Basis
    ref::BL
    out::BR
end
basis_l(b::ChoiBasis) = b.ref
basis_r(b::ChoiBasis) = b.out
Base.:(==)(b1::ChoiBasis, b2::ChoiBasis) = (b1.ref == b2.ref && b1.out == b2.out)
dimension(b::ChoiBasis) = dimension(b.ref)*dimension(b.out)
tensor(b1::ChoiBasis, b2::ChoiBasis) = ChoiBasis(tensor(b1.ref,b2.ref), tensor(b1.out, b2.out))

function _check_is_spinhalfbasis(b)
    for i=1:length(b)
        (b[i] isa SpinBasis && dimension(b[i]) == 2) || throw(ArgumentError("Must have only SpinBasis(1//2) to be compatible with pauli representation"))
    end
end
"""
    PauliBasis(N)

The standard Pauli operator basis for an `N` qubit space. This consists of
tensor products of the Pauli matrices I, X, Y, Z, in that order for each qubit.
The dimension of the basis is 2²ᴺ.
"""
struct PauliBasis{T<:Integer} <: Basis
    N::T
end
function PauliBasis(bl::Basis, br::Basis)
    bl == br || throw(ArgumentError("Both bases must be equal")) 
    _check_is_spinhalfbasis(bl)
    PauliBasis(length(bl))
end
basis_l(b::PauliBasis) = SpinBasis(1//2)^b.N
basis_r(b::PauliBasis) = SpinBasis(1//2)^b.N
Base.:(==)(b1::PauliBasis, b2::PauliBasis) = b1.N == b2.N
dimension(b::PauliBasis) = 4^b.N
tensor(b1::PauliBasis, b2::PauliBasis) = PauliBasis(b1.N+b2.N)

"""
    ChiBasis(N)

The basis for a Chi process matrix, which is just the Choi state in the Pauli
operator basis. However we do not use the `ChoiBasis`, partly to have easier
dispatch on types, and partly because there's no sensible way to distingish
between the "reference" and "output" systems as that information is lost in the
computational to Pauli basis transformation (i.e. two indices into one).

TODO explain better why dimension base is 2, see sec III.E.
"""
struct ChiBasis{T<:Integer} <: Basis
    Nl::T
    Nr::T
end
ChiBasis(N) = ChiBasis(N,N)
function ChiBasis(bl::Basis, br::Basis)
    _check_is_spinhalfbasis(bl)
    _check_is_spinhalfbasis(br)
    ChiBasis(length(bl), length(br))
end
basis_l(b::ChiBasis) = SpinBasis(1//2)^b.Nl
basis_r(b::ChiBasis) = SpinBasis(1//2)^b.Nr
Base.:(==)(b1::ChiBasis, b2::ChiBasis) = (b1.Nl == b2.Nl && b1.Nr == b2.Nr)
dimension(b::ChiBasis) = 2^(b.Nl+b.Nr)
tensor(b1::ChiBasis, b2::ChiBasis) = ChiBasis(b1.Nl+b2.Nl, b1.Nr+b2.Nr)

"""
    HWPauliBasis(N)

The Hesienberg-Weyl Pauli operator basis consisting of the N represents the
underlying Hilbert space dimension, not the operator basis dimension. For N>2,
this representes the operator basis formed by the generalized Pauli matrices,
also called the clock and shift matrices. The ordering is the usual one: when
the index is written in base-N and thus has only two digits, the least
significant bit gives powers of Z (the clock matrix), and most significant bit
gives powers of X (the shfit matrix).
"""
struct HWPauliBasis{T<:Integer} <: Basis
    shape::Vector{T}
end
HWPauliBasis(N::Integer) = HWPauliBasis([N])
Base.:(==)(b1::HWPauliBasis, b2::HWPauliBasis) = b1.shape == b2.shape
shape(b::HWPauliBasis) = [n^2 for n in b.shape]
dimension(b::HWPauliBasis) = prod([n^2 for n in b.shape])
Base.length(b::HWPauliBasis) = length(b.shape)
Base.getindex(b::HWPauliBasis, i) = HWPauliBasis([b.shape[i]...])


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
    if x.N > 1
        write(stream, "^$(x.N)")
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
    write(stream, "Pauli(N=$(x.N))")
end

function show(stream::IO, x::HWPauliBasis)
    write(stream, "Pauli($(x.shape))")
end

"""
Abstract Hilbert space class for all specialized Hilbert spaces.
"""
abstract type AbstractHilbertSpace end

"""
    space(a)

Return the Hilbert space of an object.
"""
function space end

"""
    HilbertSpace(shape)

Hilbert space of a single quantum system, whose dimension is given by `shape`.

A Hilbert space of a quantum system with finite (infinite) degrees of freedom has 
finite (infinite) dimension. See https://en.wikipedia.org/wiki/Quantum_state_space 
for mathematical details.

For example, the Hilbert space of a single spin-1/2 particle is
`HilbertSpace(2)`. The Hilbert space for a single quantum harmonic oscillator
mode is `HilbertSpace(Inf)`.
"""
struct HilbertSpace{S} <: AbstractHilbertSpace 
    shape::S
end
Base.:(==)(s1::HilbertSpace, s2::HilbertSpace) = s1.shape == s2.shape
Base.length(space::HilbertSpace) = space.shape

"""
    CompositeHilbertSpace(s1, s2...)

Hilbert space for a composite system of quantum objects with individual Hilbert spaces,
whose dimensions [s1.shape, s2.shape...] are contained in `shape` and Hilbert space objects 
[s1, s2...] are contained in `spaces`.

For example, the composite Hilbert space of a spin-1/2 particle coupled to a quantum
harmonic oscillator mode is `HilbertSpace(2) ⊗ HilbertSpace(Inf)`, which returns
a `CompositeHilbertSpace([2, Inf], [HilbertSpace(2), HilbertSpace(Inf)])` object.
"""
struct CompositeHilbertSpace{S,H} <: AbstractHilbertSpace
    shape::S
    spaces::H
end
CompositeHilbertSpace(spaces) = CompositeHilbertSpace([length(s) for s ∈ spaces], spaces)
CompositeHilbertSpace(spaces::HilbertSpace...) = CompositeHilbertSpace((spaces...,))
CompositeHilbertSpace(spaces::Vector) = CompositeHilbertSpace(spaces...)

Base.:(==)(s1::CompositeHilbertSpace, s2::CompositeHilbertSpace) = all(i == j for (i, j) in zip(s1.spaces, s2.spaces))

"""
    tensor(x::HilbertSpace, y::HilbertSpace, z::HilbertSpace...)

Create a [`CompositeHilbertSpace`](@ref) from the given Hilbert spaces.

Any given `CompositeHilbertSpace` is expanded so that the resulting `CompositeHilbertSpace` never
contains another `CompositeHilbertSpace`.
"""
tensor(s1::HilbertSpace, s2::HilbertSpace) = CompositeHilbertSpace([length(s1), length(s2)], (s1, s2))
tensor(s1::CompositeHilbertSpace, s2::CompositeHilbertSpace) = CompositeHilbertSpace(vcat(s1.shape..., s2.shape...), vcat(s1.spaces..., s2.spaces...))
function tensor(s1::CompositeHilbertSpace, s2::HilbertSpace)
    shape = vcat(s1.shape, length(s2))
    spaces = vcat(s1.spaces..., s2)
    CompositeHilbertSpace(shape, spaces)
end
function tensor(s1::HilbertSpace, s2::CompositeHilbertSpace)
    shape = vcat(length(s1), s2.shape)
    spaces = vcat(s1, s2.spaces...)
    CompositeHilbertSpace(shape, spaces)
end
tensor(spaces::T...) where {T<:AbstractHilbertSpace} = reduce(tensor, spaces)
tensor(space::T) where {T<:AbstractHilbertSpace} = space

function Base.:^(space::T, N::Integer) where {T<:AbstractHilbertSpace}
    if N < 1
        throw(ArgumentError("Power of a Hilbert space is only defined for positive integers."))
    end
    tensor([space for i=1:N]...)
end

"""
Exception that should be raised for an illegal algebraic operation between quantum objects
with different Hilbert spaces.
"""
mutable struct IncompatibleSpaces <: Exception end

"""
    samespaces(a, b)

Test if two objects have the same Hilbert spaces.
"""
samespaces(s1::AbstractHilbertSpace, s2::AbstractHilbertSpace) = s1 == s2

"""
    check_samespaces(a, b)

Throw an [`IncompatibleSpaces`](@ref) error if the objects don't have
the same Hilbert spaces.
"""
check_samespaces(s1::AbstractHilbertSpace, s2::AbstractHilbertSpace) = samespaces(s1, s2) || throw(IncompatibleSpaces())

function ptrace(space::CompositeHilbertSpace, indices)
    J = [i for i in 1:length(space.spaces) if i ∉ indices]
    length(J) > 0 || throw(ArgumentError("Tracing over all indices is not allowed in `ptrace`."))
    return reduced(space, J)
end
function reduced(space::CompositeHilbertSpace, indices)
    if length(indices)==0
        throw(ArgumentError("At least one subsystem must be specified in reduced."))
    elseif length(indices)==1
        return space.spaces[indices[1]]
    else
        return CompositeHilbertSpace(space.shape[indices], space.spaces[indices])
    end
end


##
# show methods
##

function show(io::IO, x::HilbertSpace)
    print(io, "HilbertSpace(shape=$(x.shape))")
end

function show(io::IO, x::CompositeHilbertSpace)
    print(io, "[")
    @inbounds for i in 1:length(x.spaces)
        print(io, x.spaces[i])
        if i != length(x.spaces)
            print(io, " ⊗ ")
        end
    end
    print(io, "]")
end
"""
    express(obj, repr::AbstractRepresentation)
    express(obj, repr::AbstractRepresentation, use::AbstractUse)

Translate a quantum object `obj` to a backend representation `repr`. It is relevant to define `use`
for formalism-specific cases, e.g., for `QuantumCliffordRepr`.
"""
function express end

"""An abstract type for the supported representation of quantum objects."""
abstract type AbstractRepresentation end
abstract type AbstractUse end
struct UseAsState <: AbstractUse end
struct UseAsOperation <: AbstractUse end
struct UseAsObservable <: AbstractUse end

express(obj) = express(obj, QuantumOpticsRepr()) # The default representation
express(s::Number, repr::AbstractRepresentation, use::AbstractUse) = s
express(s, repr::AbstractRepresentation) = express(s, repr, UseAsState())

##
# Commonly used representations -- interfaces for each one defined in separate packages
##

"""Representation using kets, bras, density matrices, and superoperators governed by `QuantumOptics.jl`.

When `lazy=true`, composite symbolic expressions (sums, products, tensor products)
are expressed into structure-preserving lazy operators (`LazySum`, `LazyProduct`,
`LazyTensor`) instead of eagerly-materialised dense/sparse operators."""
Base.@kwdef struct QuantumOpticsRepr <: AbstractRepresentation
    cutoff::Int = 2
    lazy::Bool = false
end
# preserve the historical positional one-argument constructor QuantumOpticsRepr(cutoff)
QuantumOpticsRepr(cutoff::Int) = QuantumOpticsRepr(cutoff, false)
"""Similar to `QuantumOpticsRepr`, but using trajectories instead of superoperators."""
struct QuantumMCRepr <: AbstractRepresentation end
"""Representation using tableaux governed by `QuantumClifford.jl`"""
struct CliffordRepr <: AbstractRepresentation end
"""Representation using Gaussian phase space formalism governed by `Gabs.jl`."""
struct GabsRepr{B} <: AbstractRepresentation 
    basis::B
end
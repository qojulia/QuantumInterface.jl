function express end

express(obj) = obj
express(obj) = express(obj, QuantumOpticsRepr()) # The default representation

"""An abstract type for the supported representation of quantum objects."""
abstract type AbstractRepresentation end

##
# Commonly used representations -- interfaces for each one defined in separate packages
##

"""Representation using kets, bras, density matrices, and superoperators governed by `QuantumOptics.jl`."""
@kwdef struct QuantumOpticsRepr <: AbstractRepresentation 
    cutoff::Int = 2
end
"""Similar to `QuantumOpticsRepr`, but using trajectories instead of superoperators."""
struct QuantumMCRepr <: AbstractRepresentation end
"""Representation using tableaux governed by `QuantumClifford.jl`"""
struct CliffordRepr <: AbstractRepresentation end
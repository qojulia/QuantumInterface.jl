module QuantumInterface

##
# Basis specific
##

"""
    basis(a)

Return the basis of an object.

If it's ambiguous, e.g. if an operator has a different left and right basis,
an [`IncompatibleBases`](@ref) error is thrown.
"""
function basis end

"""
Exception that should be raised for an illegal algebraic operation.
"""
mutable struct IncompatibleBases <: Exception end


##
# Standard methods
##

function apply! end

function dagger end

"""
    directsum(x, y, z...)

Direct sum of the given objects. Alternatively, the unicode
symbol ⊕ (\\oplus) can be used.
"""
function directsum end
const ⊕ = directsum
directsum() = GenericBasis(0)

function dm end

function embed end

function entanglement_entropy end

function expect end

function identityoperator end

function permutesystems end

function projector end

function project! end

function projectrand! end

function ptrace end

function reduced end

"""
    tensor(x, y, z...)

Tensor product of the given objects. Alternatively, the unicode
symbol ⊗ (\\otimes) can be used.
"""
function tensor end
const ⊗ = tensor
tensor() = throw(ArgumentError("Tensor function needs at least one argument."))

function tensor_pow end # TODO should Base.^ be the same as tensor_pow?

function traceout! end

function variance end

##
# Qubit specific
##

function nqubits end

function projectX! end

function projectY! end

function projectZ! end

function projectXrand! end

function projectYrand! end

function projectZrand! end

function reset_qubits! end

##
# Quantum optics specific
##

function coherentstate end

function thermalstate end

function displace end

function squeeze end

function wigner end


include("abstract_types.jl")
include("bases.jl")
include("show.jl")

include("linalg.jl")
include("tensor.jl")
include("embed_permute.jl")
include("expect_variance.jl")
include("identityoperator.jl")

include("julia_base.jl")
include("julia_linalg.jl")
include("sparse.jl")

include("sortedindices.jl")

end # module

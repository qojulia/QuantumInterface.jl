module QuantumInterface

function apply! end

function dagger end

function directsum end
const ⊕ = directsum

function dm end

function embed end

function entanglement_entropy end

function expect end

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

##
# Qubit specific
#

function nqubits end

function projectX! end

function projectY! end

function projectZ! end

function projectXrand! end

function projectYrand! end

function projectZrand! end

function reset_qubits! end

##
# Bases
##

include("bases.jl")
include("abstract_types.jl")

end # module

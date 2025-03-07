module QuantumInterface

import Base: ==, +, -, *, /, ^, length, one, exp, conj, conj!, transpose, copy
import LinearAlgebra: tr, ishermitian, norm, normalize, normalize!
import Base: show, summary
import SparseArrays: sparse, spzeros, AbstractSparseMatrix # TODO move to an extension

function apply! end

function dagger end

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

##
# Metrics
##

function entropy_vn end

function fidelity end

function logarithmic_negativity end

include("bases.jl")
include("abstract_types.jl")

include("linalg.jl")
include("tensor.jl")
include("embed_permute.jl")
include("expect_variance.jl")
include("identityoperator.jl")

include("julia_base.jl")
include("julia_linalg.jl")
include("sparse.jl")

include("sortedindices.jl")
include("express.jl")

end # module

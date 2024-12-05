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
    fullbasis(a)

Return the full basis of an object.

Returns subtype of `Basis` when a is a subtype of `StateVector`.
Returns a subtype of `OperatorBasis` a is a subtype of `AbstractOperator`.
Returns a subtype of `SuperOperatorBasis` when a is a subtype of `AbstractSuperOperator`.
"""
function fullbasis end

"""
    length(b::Basis)

Total dimension of the Hilbert space.
"""
function length end

function bases end

function spinnumber end

function cutoff end

function offset end

##
# Standard methods
##

"""
    multiplicable(a, b)

Check if any two subtypes of `StateVector`, `AbstractOperator`,
or `AbstractSuperOperator` can be multiplied in the given order.

Spcefically this checks whether the right basis of a is equal
to the left basis of b
"""
function multiplicable end

"""
    check_multiplicable(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are
not multiplicable as determined by `multiplicable(a, b)`.

If the macro `@compatiblebases` is used anywhere up the call stack,
this check is disabled.
"""
function check_multiplicable end

"""
    addible(a, b)

Check if any two subtypes of `StateVector`, `AbstractOperator`,
or `AbstractSuperOperator` can be added together.

Spcefically this checks whether the left basis of a is equal
to the left basis of b and whether the right basis of a is equal
to the right basis of b.
"""
function addible end

"""
    check_addible(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are
not addible as determined by `addible(a, b)`.

If the macro `@compatiblebases` is used anywhere up the call stack,
this check is disabled.
"""
function check_addible end

"""
    issquare(a)

Check if any two subtypes of `StateVector`, `AbstractOperator`,
or `AbstractSuperOperator` are square.

Spcefically this checks whether the left basis of a is equal
to the right basis of a.
For subtypes of `StateVector` this is always false.
"""
function addible end

"""
    check_issquare(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are
not addible as determined by `addible(a, b)`.

If the macro `@compatiblebases` is used anywhere up the call stack,
this check is disabled.
"""
function check_addible end

const BASES_CHECK = Ref(true)

"""
    @compatiblebases

Macro to skip checks for compatible bases. Useful for `*`, `expect` and similar
functions.
"""
macro compatiblebases(ex)
    return quote
        BASES_CHECK.x = false
        local val = $(esc(ex))
        BASES_CHECK.x = true
        val
    end
end

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

function equal_bases(a, b)
    Base.depwarn("`==` should be preferred over `equal_bases`!", :equal_bases)
    if a===b
        return true
    end
    for i=1:length(a)
        if a[i]!=b[i]
            return false
        end
    end
    return true
end

struct PauliBasis{S,B} <: Basis
    shape::S
    bases::B
    function PauliBasis(num_qubits::T) where {T<:Integer}
        Base.depwarn("`PauliBasis` will be removed in future versions as it does not represent the Pauli operator basis. SpinBasis(1//2)^N or NLevelBasis(2)^N should likely be used instead.", :PauliBasis)
        shape = [2 for _ in 1:num_qubits]
        bases = Tuple(SpinBasis(1//2) for _ in 1:num_qubits)
        return new{typeof(shape),typeof(bases)}(shape, bases)
    end
end

Base.:(==)(pb1::PauliBasis, pb2::PauliBasis) = length(pb1.bases) == length(pb2.bases)
Base.length(b::PauliBasis) = prod(b.shape)

# TODO: figure out how to deprecate abstract type
abstract type AbstractSuperOperator end

function basis(a::AbstractSuperOperator)
    Base.depwarn("`AbstractSuperOperator` will be removed in a future version", :equal_shape)
    (check_samebases(a); a.basis_l[1]) # FIXME issue #12
end

function check_samebases(a::AbstractSuperOperator)
    Base.depwarn("`AbstractSuperOperator` will be removed in a future version", :equal_shape)
    check_samebases(a.basis_l[1], a.basis_r[1]) # FIXME issue #12
    check_samebases(a.basis_l[2], a.basis_r[2]) # FIXME issue #12
end

function equal_shape(a, b)
    Base.depwarn("`equal_shape` will be removed in a future version", :equal_shape)
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

# TODO: figure out how to deprecate a macro
macro samebases(ex)
    return quote
        BASES_CHECK.x = false
        local val = $(esc(ex))
        BASES_CHECK.x = true
        val
    end
end

function samebases(b1::Tuple{Basis, Basis}, b2::Tuple{Basis, Basis})
    Base.depwarn("`samebases(b1:Basis, b2:Basis)`, `addible` or `multiplicable` should be used instead of `samebases(b1::Tuple{Basis, Basis}, b2::Tuple{Basis, Basis})`. ", :check_samebases)
    b1==b2 # for checking superoperators
end

function samebases(a::AbstractOperator)
    Base.depwarn("`multiplicable(a,a)` should be used instead of `samebases(a::AbstractOperator)`. ", :check_samebases)
    samebases(a.basis_l, a.basis_r)::Bool # FIXME issue #12
end

function samebases(a::AbstractOperator, b::AbstractOperator)
    Base.depwarn("``addible(a,b)` should be used instead of `samebases(a::AbstractOperator, b::AbstractOperator)`. ", :check_samebases)
    samebases(a.basis_l, b.basis_l)::Bool && samebases(a.basis_r, b.basis_r)::Bool # FIXME issue #12
end

function check_samebases(a::Union{AbstractOperator, AbstractSuperOperator})
    Base.depwarn("`check_multiplicable(a,a)` should be used instead of `check_samebases(a::Union{AbstractOperator, AbstractSuperOperator})`. ", :check_samebases)
    check_samebases(a.basis_l, a.basis_r) # FIXME issue #12
end

function multiplicable(b1::Basis, b2::Basis)
    Base.depwarn("`==` between bases should be used instead of `multiplicable`.", :check_samebases)
    b1==b2
end

function multiplicable(b1::CompositeBasis, b2::CompositeBasis)
    Base.depwarn("`==` between bases should be used instead of `multiplicable`.", :check_samebases)
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

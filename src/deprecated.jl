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
        Base.depwarn("`PauliBasis` is currently not the Pauli operator basis and its behavior will be changed in future versions.", :PauliBasis)
        shape = [2 for _ in 1:num_qubits]
        bases = Tuple(SpinBasis(1//2) for _ in 1:num_qubits)
        return new{typeof(shape),typeof(bases)}(shape, bases)
    end
end

Base.:(==)(pb1::PauliBasis, pb2::PauliBasis) = length(pb1.bases) == length(pb2.bases)


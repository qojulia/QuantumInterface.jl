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

Base.@deprecate PauliBasis(num_qubits) NQubitBasis(num_qubits) false

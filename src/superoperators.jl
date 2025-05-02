import LinearAlgebra: vec

"""
    vec(op)

Create a vectorized operator compatible with application with superoperators.
"""
vec(op::AbstractOperator) = throw(ArgumentError("vec() is not defined for this type of operator: $(typeof(op))."))

"""
    unvec(op)

Converta a vectorized operator to a normal operator.
"""
unvec(op::AbstractKet) = throw(ArgumentError("unvec() is not defined for this type of operator: $(typeof(op))."))

"""
    super(op)

Converts to superoperator representation
"""
super(op::AbstractOperator) = throw(ArgumentError("super() is not defined for this type of operator: $(typeof(op))."))

"""
    choi(op)

Converts to choi state representation
"""
choi(op::AbstractOperator) = throw(ArgumentError("choi() is not defined for this type of operator: $(typeof(op))."))

"""
    kraus(op)

Converts to kraus operator sum representation
"""
kraus(op::AbstractOperator) = throw(ArgumentError("kraus() is not defined for this type of operator: $(typeof(op))."))
kraus(op::Vector{AbstractOperator}) = throw(ArgumentError("kraus() is not defined for this type of operator: $(typeof(op))."))

"""
    stinespring(op)

Converts to (Pauli) chi process matrix representation
"""
stinespring(op::AbstractOperator) = throw(ArgumentError("stinespring() is not defined for this type of operator: $(typeof(op))."))


"""
    pauli(op)

Converts to pauli vector and transfer matrix representation
"""
pauli(op::AbstractOperator) = throw(ArgumentError("pauli() is not defined for this type of operator: $(typeof(op))."))

"""
    chi(op)

Converts to (Pauli) chi process matrix representation
"""
chi(op::AbstractOperator) = throw(ArgumentError("chi() is not defined for this type of operator: $(typeof(op))."))

"""
    spre(op)

Create a super-operator equivalent for right side operator multiplication.

For operators ``A``, ``B`` the relation

```math
    \\mathrm{spre}(A) B = A B
```

holds. `op` can be a dense or a sparse operator.
"""
function spre(op::AbstractOperator)
    multiplicable(op, op) || throw(ArgumentError("It's not clear what spre of a non-square operator should be. See issue #113"))
    sprepost(op, identityoperator(op))
end

"""
    spost(op)

Create a super-operator equivalent for left side operator multiplication.

For operators ``A``, ``B`` the relation

```math
    \\mathrm{spost}(A) B = B A
```

holds. `op` can be a dense or a sparse operator.
"""
function spost(op::AbstractOperator)
    multiplicable(op, op) || throw(ArgumentError("It's not clear what spost of a non-square operator should be. See issue #113"))
    sprepost(identityoperator(op), op)
end

"""
    sprepost(op)

Create a super-operator equivalent for left and right side operator multiplication.

For operators ``A``, ``B``, ``C`` the relation

```math
    \\mathrm{sprepost}(A, B) C = A C B
```

holds. `A` ond `B` can be dense or a sparse operators.
"""
sprepost(A::AbstractOperator, B::AbstractOperator) = throw(ArgumentError("sprepost() is not defined for these types of operator: $(typeof(A)) and  $(typeof(B))."))


function _check_input(H::AbstractOperator, J::Vector, Jdagger::Vector, rates)
    for j=J
        @assert isa(j, AbstractOperator)
    end
    for j=Jdagger
        @assert isa(j, AbstractOperator)
    end
    @assert length(J)==length(Jdagger)
    if isa(rates, Matrix{<:Number})
        @assert size(rates, 1) == size(rates, 2) == length(J)
    elseif isa(rates, Vector{<:Number})
        @assert length(rates) == length(J)
    end
end

"""
    liouvillian(H, J; rates, Jdagger)

Create a super-operator equivalent to the master equation so that ``\\dot ρ = S ρ``.

The super-operator ``S`` is defined by

```math
S ρ = -\\frac{i}{ħ} [H, ρ] + \\sum_i J_i ρ J_i^† - \\frac{1}{2} J_i^† J_i ρ - \\frac{1}{2} ρ J_i^† J_i
```

# Arguments
* `H`: Hamiltonian.
* `J`: Vector containing the jump operators.
* `rates`: Vector or matrix specifying the coefficients for the jump operators.
* `Jdagger`: Vector containing the hermitian conjugates of the jump operators. If they
             are not given they are calculated automatically.
"""
function liouvillian(H, J; rates=ones(length(J)), Jdagger=dagger.(J))
    _check_input(H, J, Jdagger, rates)
    L = spre(-1im*H) + spost(1im*H)
    if isa(rates, AbstractMatrix)
        for i=1:length(J), j=1:length(J)
            jdagger_j = rates[i,j]/2*Jdagger[j]*J[i]
            L -= spre(jdagger_j) + spost(jdagger_j)
            L += spre(rates[i,j]*J[i]) * spost(Jdagger[j])
        end
    elseif isa(rates, AbstractVector)
        for i=1:length(J)
            jdagger_j = rates[i]/2*Jdagger[i]*J[i]
            L -= spre(jdagger_j) + spost(jdagger_j)
            L += spre(rates[i]*J[i]) * spost(Jdagger[i])
        end
    end
    return L
end

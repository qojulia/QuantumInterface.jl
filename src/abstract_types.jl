"""
Abstract type for `Bra` and `Ket` states.

The state vector type stores an abstract state with respect to a certain
Hilbert space basis.
All deriving types must define the `fullbasis` function which
returns the state vector's underlying `Basis`.
"""
abstract type StateVector{B<:Basis} end
abstract type AbstractKet{B} <: StateVector{B} end
abstract type AbstractBra{B} <: StateVector{B} end

"""
Abstract type for all operators which represent linear maps between two
Hilbert spaces with respect to a given basis in each space.

All deriving operator types must define the `fullbasis` function which
returns the operator's underlying `OperatorBasis`.

For fast time evolution also at least the function
`mul!(result::Ket,op::AbstractOperator,x::Ket,alpha,beta)` should be
implemented. Many other generic multiplication functions can be defined in
terms of this function and are provided automatically.

See [TODO: reference operators.md in docs]
"""
abstract type AbstractOperator{B<:OperatorBasis} end

"""
Abstract type for all super-operators which represent linear maps between two
operator spaces with respect to a given basis for each space.

All deriving operator types must define the `fullbasis` function which
returns the operator's underlying `SuperOperatorBasis`.

See [TODO: reference superoperators.md in docs]
```
"""
abstract type AbstractSuperOperator{B<:SuperOperatorBasis} end

function summary(stream::IO, x::AbstractOperator)
    b = fullbasis(x)
    print(stream, "$(typeof(x).name.name)(dim=$(length(b.left))x$(length(b.right)))\n")
    if samebases(b)
        print(stream, "  basis: ")
        show(stream, basis(b))
    else
        print(stream, "  basis left:  ")
        show(stream, b.left)
        print(stream, "\n  basis right: ")
        show(stream, b.right)
    end
end

show(stream::IO, x::AbstractOperator) = summary(stream, x)

traceout!(s::StateVector, i) = ptrace(s,i)

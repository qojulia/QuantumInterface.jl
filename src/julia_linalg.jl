"""
    ishermitian(op::AbstractOperator)

Check if an operator is Hermitian.
"""
ishermitian(op::AbstractOperator) = arithmetic_unary_error(ishermitian, op)

"""
    tr(x::AbstractOperator)

Trace of the given operator.
"""
tr(x::AbstractOperator) = arithmetic_unary_error("Trace", x)

"""
    norm(x::StateVector)

Norm of the given bra or ket state.
"""
norm(x::StateVector) = norm(x.data)

"""
    normalize(x::StateVector)

Return the normalized state so that `norm(x)` is one.
"""
normalize(x::StateVector) = x/norm(x)

"""
    normalize!(x::StateVector)

In-place normalization of the given bra or ket so that `norm(x)` is one.
"""
normalize!(x::StateVector) = (normalize!(x.data); x)

"""
    normalize(op)

Return the normalized operator so that its `tr(op)` is one.
"""
normalize(op::AbstractOperator) = op/tr(op)

"""
    normalize!(op)

In-place normalization of the given operator so that its `tr(x)` is one.
"""
normalize!(op::AbstractOperator) = throw(ArgumentError("normalize! is not defined for this type of operator: $(typeof(op)).\n You may have to fall back to the non-inplace version 'normalize()'."))

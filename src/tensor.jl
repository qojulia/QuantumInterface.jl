"""
    tensor(x::AbstractOperator, y::AbstractOperator, z::AbstractOperator...)

Tensor product ``\\hat{x}⊗\\hat{y}⊗\\hat{z}⊗…`` of the given operators.
"""
tensor(a::AbstractOperator, b::AbstractOperator) = arithmetic_binary_error("Tensor product", a, b)
tensor(op::AbstractOperator) = op
tensor(operators::AbstractOperator...) = reduce(tensor, operators)
tensor(state::StateVector) = state
tensor(states::Vector{T}) where T<:StateVector = reduce(tensor, states)

"""
    tensor_pow(a, N)

Gives the tensor product of `a` `N` times.
"""
tensor_pow(a, N) = Base.power_by_squaring(a, N; mul=⊗)

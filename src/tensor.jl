"""
    tensor(x::AbstractOperator, y::AbstractOperator, z::AbstractOperator...)

Tensor product ``\\hat{x}⊗\\hat{y}⊗\\hat{z}⊗…`` of the given operators.
"""
tensor(a::AbstractOperator, b::AbstractOperator) = arithmetic_binary_error("Tensor product", a, b)
tensor(op::AbstractOperator) = op
tensor(operators::AbstractOperator...) = reduce(tensor, operators)
tensor(state::StateVector) = state
tensor(states::StateVector...) = reduce(tensor, states)
tensor(states::Vector{T}) where T<:StateVector = reduce(tensor, states)

"""
    tensor_pow(a, N)

Gives the tensor product of `a` `N` times.
"""
tensor_pow(a, N) = tensor_pow_by_squaring(a, N)

# Copied from Base.power_by_squaring as `mul` keyword dosn't work without implementing `*`
function tensor_pow_by_squaring(x, p::Integer)
    if p == 1
        return x
    elseif p == 2
        return tensor(x, x)
    elseif p < 1
        throw(DomainError("Cannot take tensor_pow to power less than one"))
    end
    t = trailing_zeros(p) + 1
    p >>= t
    while (t -= 1) > 0
        x = tensor(x, x)
    end
    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) >= 0
            x = tensor(x, x)
        end
        y = tensor(y, x)
    end
    return y
end

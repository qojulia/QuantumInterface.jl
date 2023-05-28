"""
    expect(index, op, state)

If an `index` is given, it assumes that `op` is defined in the subsystem specified by this number.
"""
function expect(indices, op::AbstractOperator{B1,B2}, state::AbstractOperator{B3,B3}) where {B1,B2,B3<:CompositeBasis}
    N = length(state.basis_l.shape)
    indices_ = complement(N, indices)
    expect(op, ptrace(state, indices_))
end

expect(index::Integer, op::AbstractOperator{B1,B2}, state::AbstractOperator{B3,B3}) where {B1,B2,B3<:CompositeBasis} = expect([index], op, state)
expect(op::AbstractOperator, states::Vector) = [expect(op, state) for state=states]
expect(indices, op::AbstractOperator, states::Vector) = [expect(indices, op, state) for state=states]

expect(op::AbstractOperator{B1,B2}, state::AbstractOperator{B2,B2}) where {B1,B2} = tr(op*state)

"""
    variance(index, op, state)

If an `index` is given, it assumes that `op` is defined in the subsystem specified by this number
"""
function variance(indices, op::AbstractOperator{B,B}, state::AbstractOperator{BC,BC}) where {B,BC<:CompositeBasis}
    N = length(state.basis_l.shape)
    indices_ = complement(N, indices)
    variance(op, ptrace(state, indices_))
end

variance(index::Integer, op::AbstractOperator{B,B}, state::AbstractOperator{BC,BC}) where {B,BC<:CompositeBasis} = variance([index], op, state)
variance(op::AbstractOperator, states::Vector) = [variance(op, state) for state=states]
variance(indices, op::AbstractOperator, states::Vector) = [variance(indices, op, state) for state=states]

function variance(op::AbstractOperator{B,B}, state::AbstractOperator{B,B}) where B
    expect(op*op, state) - expect(op, state)^2
end

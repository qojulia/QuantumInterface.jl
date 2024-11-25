samebases(a::AbstractOperator) = samebases(fullbasis(a))::Bool
samebases(a::AbstractSuperOperator) = samebases(fullbasis(a))::Bool
samebases(a::AbstractOperator, b::AbstractOperator) = samebases(fullbasis(a), fullbasis(b))::Bool
samebases(a::AbstractSuperOperator, b::AbstractSuperOperator) = samebases(fullbasis(a), fullbasis(b))::Bool
check_samebases(a::AbstractOperator) = check_samebases(fullbasis(a))::Bool
check_samebases(a::AbstractSuperOperator) = check_samebases(fullbasis(a))::Bool
check_samebases(a::AbstractOperator, b::AbstractOperator) = check_samebases(fullbasis(a), fullbasis(b))::Bool
check_samebases(a::AbstractSuperOperator, b::AbstractSuperOperator) = check_samebases(fullbasis(a), fullbasis(b))::Bool

multiplicable(a::AbstractOperator, b::AbstractOperator) = multiplicable(fullbasis(a).right, fullbasis(b).left)
dagger(a::AbstractOperator) = arithmetic_unary_error("Hermitian conjugate", a)
transpose(a::AbstractOperator) = arithmetic_unary_error("Transpose", a)
directsum(a::AbstractOperator...) = reduce(directsum, a)
ptrace(a::AbstractOperator, index) = arithmetic_unary_error("Partial trace", a)
_index_complement(b::CompositeBasis, indices) = complement(length(b.bases), indices)
reduced(a, indices) = ptrace(a, _index_complement(basis(a), indices))

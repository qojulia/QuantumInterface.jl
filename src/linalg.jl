samebases(a::AbstractOperator) = samebases(a.basis_l, a.basis_r)::Bool # FIXME issue #12
samebases(a::AbstractOperator, b::AbstractOperator) = samebases(a.basis_l, b.basis_l)::Bool && samebases(a.basis_r, b.basis_r)::Bool # FIXME issue #12
check_samebases(a::Union{AbstractOperator, AbstractSuperOperator}) = check_samebases(a.basis_l, a.basis_r) # FIXME issue #12
multiplicable(a::AbstractOperator, b::AbstractOperator) = multiplicable(a.basis_r, b.basis_l) # FIXME issue #12
directsum(a::AbstractOperator...) = reduce(directsum, a)
ptrace(a::AbstractOperator, index) = arithmetic_unary_error("Partial trace", a)
_index_complement(b::CompositeBasis, indices) = complement(length(b.bases), indices)
reduced(a, indices) = ptrace(a, _index_complement(basis(a), indices))

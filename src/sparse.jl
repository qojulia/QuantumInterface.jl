# TODO make an extension?

# dense(a::AbstractOperator) = arithmetic_unary_error("Conversion to dense", a)

"""
    sparse(op::AbstractOperator)

Convert an arbitrary operator into a [`SparseOperator`](@ref).
"""
sparse(a::AbstractOperator) = throw(ArgumentError("Direct conversion from $(typeof(a)) not implemented. Use sparse(full(op)) from QuantumOptics instead."))

function ptrace(x::AbstractSparseMatrix, shape_nd, indices) # A fairly QuantumOptics-specific method, left here to avoid piracy
    shape_nd = (shape_nd...,)
    N = div(length(shape_nd), 2)
    shape_2d = (x.m, x.n)
    shape_nd_after = ([i ∈ indices || i-N ∈ indices ? 1 : shape_nd[i] for i=1:2*N]...,)
    shape_2d_after = (prod(shape_nd_after[1:N]), prod(shape_nd_after[N+1:end]))
    I_nd_after_max = CartesianIndex(shape_nd_after...)
    y = spzeros(eltype(x), shape_2d_after...)
    for I in eachindex(x)::CartesianIndices{2} # Manual type assertions to help JET
        println(I.I)
        I_nd = sub2sub(shape_2d, shape_nd, I)::CartesianIndex{2}
        if I_nd.I[indices] != I_nd.I[indices .+ N]
            continue
        end
        I_after = sub2sub(shape_nd_after, shape_2d_after, min(I_nd, I_nd_after_max)::CartesianIndex{2})
        y[I_after] += x[I]
    end
    y
end

function sub2sub(shape1, shape2, I::CartesianIndex{2})::CartesianIndex{2} # Manual type assertions to help JET
    linearindex = LinearIndices(shape1)[I.I...]
    CartesianIndices(shape2)[linearindex]
end

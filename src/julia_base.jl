# Common error messages
# We use these error messages instead of:
# - MethodErrors, which are not very informative about the fact that there might very well be a convenient conversion to the desired type
# - error hints, which are a bit too cumbersome to register for so many different types and functions
# TODO use error hints instead of these custom error messages
arithmetic_unary_error(funcname, x::AbstractOperator) = throw(ArgumentError("$funcname is not defined for this type of operator: $(typeof(x)).\nYou can try to convert to another operator type first, e.g. with `dense()` or `sparse()` from `QuantumOptics` or other type conversions from other packages."))
arithmetic_binary_error(funcname, a::AbstractOperator, b::AbstractOperator) = throw(ArgumentError("$funcname is not defined for this combination of types of operators: $(typeof(a)), $(typeof(b)).\nYou can try to convert to another (common) operator type first, e.g. with `dense()` or `sparse()` from `QuantumOptics` or other type conversions from other packages."))
addnumbererror() = throw(ArgumentError("Can't add or subtract a number and an operator. You probably want `op + identityoperator(op)*x`."))


##
# States
##

-(a::StateVector) = arithmetic_unary_error("Negation", a)
*(a::StateVector, b::Number) = b*a
length(a::StateVector) = length(a.basis)::Int # FIXME issue #12
basis(a::StateVector) = a.basis # FIXME issue #12
directsum(x::StateVector...) = reduce(directsum, x)

# Array-like functions
Base.size(x::StateVector) = length(basis(x))
Base.ndims(x::StateVector) = 1
Base.ndims(::Type{<:StateVector}) = 1

# Broadcasting
Base.broadcastable(x::StateVector) = x

Base.adjoint(a::StateVector) = dagger(a)


##
# Operators
##

length(a::AbstractOperator) = length(a.basis_l)::Int*length(a.basis_r)::Int  # FIXME issue #12
basis(a::AbstractOperator) = (check_samebases(a); a.basis_l) # FIXME issue #12
basis(a::AbstractSuperOperator) = (check_samebases(a); a.basis_l[1]) # FIXME issue #12

# Ensure scalar broadcasting
Base.broadcastable(x::AbstractOperator) = Ref(x)

# Arithmetic operations
+(a::AbstractOperator, b::AbstractOperator) = arithmetic_binary_error("Addition", a, b)
+(a::Number, b::AbstractOperator) = addnumbererror()
+(a::AbstractOperator, b::Number) = addnumbererror()
+(a::AbstractOperator) = a

-(a::AbstractOperator) = arithmetic_unary_error("Negation", a)
-(a::AbstractOperator, b::AbstractOperator) = arithmetic_binary_error("Subtraction", a, b)
-(a::Number, b::AbstractOperator) = addnumbererror()
-(a::AbstractOperator, b::Number) = addnumbererror()

*(a::AbstractOperator, b::AbstractOperator) = arithmetic_binary_error("Multiplication", a, b)
^(a::AbstractOperator, b::Integer) = Base.power_by_squaring(a, b)

"""
    exp(op::AbstractOperator)

Operator exponential.
"""
exp(op::AbstractOperator) = throw(ArgumentError("exp() is not defined for this type of operator: $(typeof(op)).\nTry to convert to dense operator first with dense()."))

Base.size(op::AbstractOperator) = (length(op.basis_l),length(op.basis_r)) # FIXME issue #12
function Base.size(op::AbstractOperator, i::Int)
    i < 1 && throw(ErrorException("dimension index is < 1"))
    i > 2 && return 1
    i==1 ? length(op.basis_l) : length(op.basis_r) # FIXME issue #12
end

Base.adjoint(a::AbstractOperator) = dagger(a)

conj(a::AbstractOperator) = arithmetic_unary_error("Complex conjugate", a)
conj!(a::AbstractOperator) = conj(a::AbstractOperator)

# Common error messages
arithmetic_unary_error(funcname, x::AbstractOperator) = throw(ArgumentError("$funcname is not defined for this type of operator: $(typeof(x)).\nTry to convert to another operator type first with e.g. dense() or sparse()."))
arithmetic_binary_error(funcname, a::AbstractOperator, b::AbstractOperator) = throw(ArgumentError("$funcname is not defined for this combination of types of operators: $(typeof(a)), $(typeof(b)).\nTry to convert to a common operator type first with e.g. dense() or sparse()."))
addnumbererror() = throw(ArgumentError("Can't add or subtract a number and an operator. You probably want 'op + identityoperator(op)*x'."))


##
# States
##

-(a::T) where {T<:StateVector} = T(a.basis, -a.data) # FIXME
*(a::StateVector, b::Number) = b*a
copy(a::T) where {T<:StateVector} = T(a.basis, copy(a.data)) # FIXME
length(a::StateVector) = length(basis(a))::Int
basis(a::StateVector) = fullbasis(a)
directsum(x::StateVector...) = reduce(directsum, x)

# Array-like functions
Base.size(x::StateVector) = size(x.data) # FIXME
@inline Base.axes(x::StateVector) = axes(x.data) #FIXME
Base.ndims(x::StateVector) = 1
Base.ndims(::Type{<:StateVector}) = 1
Base.eltype(x::StateVector) = eltype(x.data) # FIXME

# Broadcasting
Base.broadcastable(x::StateVector) = x

Base.adjoint(a::StateVector) = dagger(a)


##
# Operators
##

length(a::AbstractOperator) = (b=fullbasis(a); length(b.left)::Int*length(b.right)::Int)
basis(a::AbstractOperator) = (b=fullbasis(a); check_samebases(b); b.left)
basis(a::AbstractSuperOperator) = (b=fullbasis(a); check_samebases(b); b.left.left)

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

Base.size(op::AbstractOperator) = (b=fullbasis(op); (length(b.left),length(b.right)))
function Base.size(op::AbstractOperator, i::Int)
    i < 1 && throw(ErrorException("dimension index is < 1"))
    i > 2 && return 1
    i==1 ? length(fullbasis(op).left) : length(fullbasis(op).right)
end

Base.adjoint(a::AbstractOperator) = dagger(a)

conj(a::AbstractOperator) = arithmetic_unary_error("Complex conjugate", a)
conj!(a::AbstractOperator) = conj(a::AbstractOperator)

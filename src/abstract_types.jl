"""
Abstract base class for all specialized bases.

The Basis class is meant to specify a basis of the Hilbert space of the
studied system. Besides basis specific information all subclasses must
implement a shape variable which indicates the dimension of the used
Hilbert space. For a spin-1/2 Hilbert space this would be the
vector `[2]`. A system composed of two spins would then have a
shape vector `[2 2]`.

Composite systems can be defined with help of the [`CompositeBasis`](@ref)
class.
"""
abstract type Basis end

"""
Abstract base class for `Bra` and `Ket` states.

The state vector class stores the coefficients of an abstract state
in respect to a certain basis. These coefficients are stored in the
`data` field and the basis is defined in the `basis`
field.
"""
abstract type StateVector{B,T} end
abstract type AbstractKet{B,T} <: StateVector{B,T} end
abstract type AbstractBra{B,T} <: StateVector{B,T} end

"""
Abstract base class for all operators.

All deriving operator classes have to define the fields
`basis_l` and `basis_r` defining the left and right side bases.

For fast time evolution also at least the function
`mul!(result::Ket,op::AbstractOperator,x::Ket,alpha,beta)` should be
implemented. Many other generic multiplication functions can be defined in
terms of this function and are provided automatically.
"""
abstract type AbstractOperator{BL,BR} end

"""
Base class for all super operator classes.

Super operators are bijective mappings from operators given in one specific
basis to operators, possibly given in respect to another, different basis.
To embed super operators in an algebraic framework they are defined with a
left hand basis `basis_l` and a right hand basis `basis_r` where each of
them again consists of a left and right hand basis.
```math
A_{bl_1,bl_2} = S_{(bl_1,bl_2) ↔ (br_1,br_2)} B_{br_1,br_2}
\\\\
A_{br_1,br_2} = B_{bl_1,bl_2} S_{(bl_1,bl_2) ↔ (br_1,br_2)}
```
"""
abstract type AbstractSuperOperator{B1,B2} end

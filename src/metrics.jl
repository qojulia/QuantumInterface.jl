"""
Calculate the Von Neumann entropy of a density operator, defined as

```math
S(\\rho) = -Tr(\\rho \\log(\\rho))
```

wherein it is understood that ``0 \\log(0) \\equiv 0``.

Consult specific implementation for function arguments and logarithmic basis.
"""
function entropy_vn end
# avoids causing a breaking change in QuantumOpticsBase.jl
entropy_vn(::StateVector; kwargs...) = 0
"""
Calculate the joint fidelity of two density operators, defined as

```math
F(\\rho, \\sigma) = Tr(\\sqrt{\\sqrt{\\rho} \\sigma \\sqrt{\\rho}}).
```

Consult specific implementation for function arguments.
"""
function fidelity end

"""
Calculate the logarithmic negativity of a density operator partition, defined as

```math
N(\\rho) = \\log\\|\\rho^{T_B}\\|_1
```

wherein ``\\rho^{T_B}`` denotes partial transposition as per the partition and
``\\|\\bullet\\|_1`` is the trace norm.


Consult specific implementation for function arguments and logarithmic basis.
"""
function logarithmic_negativity end

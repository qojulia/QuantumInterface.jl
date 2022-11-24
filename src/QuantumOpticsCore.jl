module QuantumOpticsCore

"""
    tensor

Generic tensor product.
"""
function tensor end
const ⊗ = tensor


"""
    create

Bosonic creation operator.
"""
function create end

"""
    destroy

Bosonic annihilation operator.
"""
function destroy end

"""
    transition

Discrete level transitions.
"""
function transition end


end  # module
using Test
using QuantumInterface: QuantumOpticsRepr

@testset "QuantumOpticsRepr constructors" begin
    # the default stays eager
    @test QuantumOpticsRepr() == QuantumOpticsRepr(2, false)
    @test QuantumOpticsRepr().cutoff == 2
    @test QuantumOpticsRepr().lazy == false

    # the positional `cutoff` constructor is preserved
    @test QuantumOpticsRepr(4) == QuantumOpticsRepr(cutoff=4)
    @test QuantumOpticsRepr(4).cutoff == 4
    @test QuantumOpticsRepr(4).lazy == false

    # the new `lazy` option is opt-in via keyword
    @test QuantumOpticsRepr(lazy=true).lazy == true
    @test QuantumOpticsRepr(cutoff=4, lazy=true) == QuantumOpticsRepr(4, true)
end

@testset "QuantumOpticsRepr cache-key equality" begin
    # `express` caches conversions in a `Dict` keyed on the representation value, so an
    # eager and a lazy representation of the same object must compare unequal (and hash
    # differently) or the cache would hand back the wrong type.
    @test QuantumOpticsRepr(4) == QuantumOpticsRepr(cutoff=4, lazy=false)
    @test QuantumOpticsRepr(4) != QuantumOpticsRepr(4, true)
    @test hash(QuantumOpticsRepr(4)) == hash(QuantumOpticsRepr(cutoff=4, lazy=false))
    @test hash(QuantumOpticsRepr(2)) != hash(QuantumOpticsRepr(2, true))
end

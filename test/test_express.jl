using Test
using QuantumInterface: QuantumOpticsRepr

@testset "express representations" begin
    repr = QuantumOpticsRepr()
    @test repr.cutoff == 2
    @test repr.lazy == false

    repr = QuantumOpticsRepr(4)
    @test repr.cutoff == 4
    @test repr.lazy == false

    repr = QuantumOpticsRepr(cutoff=5, lazy=true)
    @test repr.cutoff == 5
    @test repr.lazy == true
end

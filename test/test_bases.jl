using Test
using QuantumInterface: tensor, ⊗, ptrace, reduced, permutesystems, multiplicable
using QuantumInterface: GenericBasis, CompositeBasis, NLevelBasis, FockBasis

@testset "basis" begin

d1 = 5
d2 = 2
d3 = 6

b1 = GenericBasis(d1)
b2 = GenericBasis(d2)
b3 = GenericBasis(d3)

@test length(b1) == d1
@test length(b2) == d2
@test b1 != b2
@test b1 != FockBasis(2)
@test b1 == b1

@test tensor(b1) == b1
comp_b1 = tensor(b1, b2)
comp_uni = b1 ⊗ b2
comp_b2 = tensor(b1, b1, b2)
@test size(comp_b1) == [d1, d2]
@test size(comp_uni) == [d1, d2]
@test size(comp_b2) == [d1, d1, d2]

@test b1^3 == CompositeBasis(b1, b1, b1)
@test (b1⊗b2)^2 == CompositeBasis(b1, b2, b1, b2)
@test_throws DomainError b1^(0)

comp_b1_b2 = tensor(comp_b1, comp_b2)
@test size(comp_b1_b2) == [d1, d2, d1, d1, d2]
@test comp_b1_b2 == CompositeBasis(b1, b2, b1, b1, b2)

@test_throws ArgumentError tensor()
@test size(comp_b2) == size(tensor(b1, comp_b1))
@test comp_b2 == tensor(b1, comp_b1)

@test_throws ArgumentError ptrace(comp_b1, [1, 2])
@test ptrace(comp_b2, [1]) == ptrace(comp_b2, [2]) == comp_b1 == ptrace(comp_b2, 1)
@test ptrace(comp_b2, [1, 2]) == ptrace(comp_b1, [1])
@test ptrace(comp_b2, [2, 3]) == ptrace(comp_b1, [2])
@test ptrace(comp_b2, [2, 3]) == reduced(comp_b2, [1])
@test_throws ArgumentError reduced(comp_b1, [])

comp1 = tensor(b1, b2, b3)
comp2 = tensor(b2, b1, b3)
@test permutesystems(comp1, [2,1,3]) == comp2

@test [b1, b2] != [b1, b3]
@test comp1 != b1 ⊗ b2 ⊗ NLevelBasis(prod(size(b3)))

end # testset

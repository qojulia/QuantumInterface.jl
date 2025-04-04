using Test
using QuantumInterface
import QuantumInterface: ⊗

s = QuantumInterface

hs2 = s.HilbertSpace(2)
hsinf = s.HilbertSpace(Inf)
tp = s.CompositeHilbertSpace([2, Inf], [s.HilbertSpace(2), s.HilbertSpace(Inf)])

@test length(hs2) == 2 && length(hsinf) == Inf
@test hs2 ⊗ hsinf == tp
@test s.CompositeHilbertSpace(hs2, hsinf) == tp
@test s.CompositeHilbertSpace((hs2, hsinf)) == tp
@test s.CompositeHilbertSpace([hs2, hsinf]) == tp

@test tp ⊗ hs2 == s.CompositeHilbertSpace([2, Inf, 2], [hs2, hsinf, hs2])
@test hs2 ⊗ tp == s.CompositeHilbertSpace([2, 2, Inf], [hs2, hs2, hsinf])
@test tp ⊗ tp == s.CompositeHilbertSpace([2, Inf, 2, Inf], [hs2, hsinf, hs2, hsinf])
@test tp ⊗ hs2 ⊗ tp == s.CompositeHilbertSpace([2, Inf, 2, 2, Inf], [hs2, hsinf, hs2, hs2, hsinf])
@test s.tensor(tp) == tp

@test_throws ArgumentError hs2^0
@test hs2^3 == hs2 ⊗ hs2 ⊗ hs2

@test s.check_samespaces(hs2, hs2) == true
@test_throws s.IncompatibleSpaces s.check_samespaces(hs2, hsinf)

@test s.ptrace(tp, 2) == hs2
@test_throws ArgumentError s.ptrace(tp, [1, 2])
@test_throws ArgumentError s.reduced(tp, [])
@test s.ptrace(hs2 ⊗ hsinf ⊗ hsinf ⊗ hs2, [2, 4]) == hs2 ⊗ hsinf
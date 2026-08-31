using Aqua
using QuantumInterface
# `tensor` is a shared interface that QuantumInterface intentionally extends.
Aqua.test_all(QuantumInterface; piracies = (; treat_as_own = [QuantumInterface.tensor]))

# News

## v0.4.1 - 2025-06-06

- Define `tensor(states::StateVector...) = reduce(tensor, states)` method.

## v0.4.0 - 2025-06-04

- **(breaking)** Move methods which access a `.data` field but are defined only on an abstract types to QuantumOpticsBase
- Declare `mutual_information`, without any implemented methods.
- Move `PositionBasis` and `MomentumBasis` here from `QuantumOpticsBase`.

## v0.3.10 - 2025-04-21

- Deprecate `PauliBasis` as it is identical to `SpinBasis(1//2)^N` but without the compatibility with `CompositeBasis`.

## v0.3.9 - 2025-04-20

- Implement general `tensor_pow` via `power_by_squaring`
- Improvements to implementation of `entropy_vn`.
- Deprecate `equal_bases`.
- Add `GabsRepr` type for representations in Gabs.

## v0.3.8 - 2025-03-08

- Introduce `metrics.jl` file for metric declarations.

## v0.3.7 - 2025-01-16

- Migrate `express` functionalities and representation types from QuantumSymbolics.

## v0.3.6 - 2024-09-08

- Add `coherentstate`, `thermalstate`, `displace`, `squeeze`, `wigner`, previously from QuantumOptics.

## v0.3.5

- Fix piracies and ambiguities in `nsubsystems` accumulated downstream in QuantumSavory.

## v0.3.4

- Documentation build fix.

## v0.3.3

- Add `nsubsystems` for computing the number of subsystems in a state.

## v0.3.2

- Cleanup - removing forgotten debug print statements.

## v0.3.1

- `identitysuperoperator` declared.

## v0.3.0

- Redo some of the `identityoperator` methods.

## v0.2.2

- Move `Base.adjoint` from `QuantumOpticsBase` to `QuantumInterface`.

## v0.2.1

- Implement `basis` method for superoperators.

## v0.2.0

- Moving much more of `QuantumOptictsBase` to `QuantumInterface`, avoiding piracy. Now `QuantumInterface` takes care of abstract state/operator types and concrete bases, while `QuantumOpticsBase` has concrete Schroedinger-style implementations of concrete state/operator types. `QuantumClifford` gives tableax-style implementations of concrete state/operator types.

## v0.1.0

- first release, bringing over interfaces from `QuantumOpticsBase`, `QuantumClifford`, `QuantumSavory`, and `QuantumSymbolics`

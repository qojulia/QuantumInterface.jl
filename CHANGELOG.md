# News

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


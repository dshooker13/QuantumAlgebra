# Changelog

## 1.0.1
- Fixed failing `Matrix-spin-1-commutator` test caused by a minor issue in `iSpinLadder` where `Conjugate[\[HBar]]!=\[HBar]`.
- Fixed minor issue in `iZeroQ`, which was causing `QAExpandDualExponential` to not truncate as expectd and fail its unit test.

## 1.0.0

First public release.

### Fixed (relative to the original script)
- **NCAlgebra function lookup.** `BeginPackage["QuantumAlgebra`"]` removed NCAlgebra's contexts from
  `$ContextPath`, so `NCE` and `NCRR` could be read as new, undefined private symbols. NCAlgebra is now
  loaded before the package starts, and its contexts are added to the package's private section.
- **NCAlgebra loading check.** It looked for `"NCalgebra`"` (wrong capitalization), so the check never matched.
- **ħ leaked as `QuantumAlgebra`Private`ħ`.** Now uses `Global`ħ` by default, with an `"HBar"` option.
- **Infinite loops.** Non-terminating series and cyclic rules used to hang forever. They now stop at a
  limit and show a warning.
- **Fragile termination tests.** `==` on symbolic expressions can stay unevaluated. Tests now use
  `SameQ`, plus `Simplify` for zero checks.
- **Powers.** Rules like `a**c -> c**a + 1` now also apply inside `a^2**c`, using the new
  word-rewriting engine.
- **`QASymmetrize`** checked `NCSymmetricQ`, which tests self-adjointness, not symmetry under swapping
  operators, so it often gave a wrong warning. That check is removed. Reordering now works on whole words.
- **`QACGTables`** printed its output and returned `Null`. It now returns the tables.
- **Examples** now use `c` instead of `ad` for the creation operator, since NCAlgebra only treats single
  letters as noncommutative by default.

### Added
- Ready-made algebras: `QABosonRules`, `QAFermionRules`, `QACanonicalRules`, `QASpinRules`,
  `QASpinLadderRules`, `QAPauliRules`, `QACommutingRules` (multi-mode supported).
- `QACheckRules`: critical-pair check for conflicting, cyclic or non-confluent rules.
- `QADagger` and `QAHermitianQ`.
- `QAHeisenbergEOM` and `QAHeisenbergSolve`, which solves exactly by linear closure and `MatrixExp`.
- Closed-form summation of repeating commutator series, with linear closure as a fallback
  (`"Expansion"` option).
- `QABCH`, `QAOperatorSubstitute`, `QADisplace`, `QASqueeze`.
- `QAFockExpectation`, `QAFockMatrixElement`, `QACoherentExpectation`, `QAThermalExpectation`.
- `QAToMatrix`, `QABosonMatrices`, `QASpinMatrices`, `QASpinLadderMatrices`, `QAPauliMatrices`.
- `QAAntiCommutator`, `QANestedCommutator`, `QAWeylOrder`, `QACGTable`, `QACoupledState`.
- `QASetDaggerFormat` / `QAClearDaggerFormat`.
- A warning when a scalar parameter (e.g. a bare `t`) is noncommutative.
- Optional `t0` argument for `QAPropagateOperator`.
- Paclet layout, test suite, examples, README, MIT license.

### Changed
- `QAPropogateOperator` was renamed to `QAPropagateOperator`. The old name still works but shows a
  deprecation warning.
- Public symbols are `Protected`.

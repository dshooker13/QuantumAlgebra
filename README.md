# QuantumAlgebra

Operator algebra for quantum mechanics in Mathematica, built on
[NCAlgebra](https://github.com/NCAlgebra/NC). With it you can:

- **Normal order** with your own rules, or with ready-made bosonic, fermionic, canonical, spin and Pauli algebras
- **Commutators**: plain, anti and nested commutators, plus a consistency check for rule sets
- **Dynamics**: Heisenberg equations of motion and exact Heisenberg-picture operators
- **Transformations**: `e^A B e^-A`, Baker–Campbell–Hausdorff, displacement, squeezing and general operator substitutions
- **Adjoints**: Hermitian conjugates and Hermiticity checks
- **Expectation values** in Fock, coherent and thermal states
- **Matrix representations** for checking results numerically, including tensor products
- **Ordering prescriptions**: symmetrization and Weyl ordering
- **Angular momentum coupling**: Clebsch–Gordan tables and coupled states

## Requirements

- Mathematica / Wolfram Engine 12.1+
- NCAlgebra 5.0+

## Installation

Clone the repo, then in Mathematica:

```wolfram
PacletDirectoryLoad["/path/to/QuantumAlgebra"];
Needs["QuantumAlgebra`"]
```

Or install it permanently:

```wolfram
PacletInstall[CreatePacletArchive["/path/to/QuantumAlgebra"]];
Needs["QuantumAlgebra`"]
```

NCAlgebra is loaded automatically if it isn't already.

## Two things to know first

**Every single lowercase letter is an operator by default.** NCAlgebra treats `a`–`z` as
noncommutative. That's convenient for operators, but scalars such as `t`, `m`, `n` or `r` must be declared:

```wolfram
SetCommutative[t, m, n, r];
```

Greek letters, `ħ` and multi-letter names are commutative unless you say otherwise.
QuantumAlgebra warns you when a parameter it expects to be a scalar (such as a time) is noncommutative.

**Create rule sets before writing expressions with multi-letter operators.** Functions like
`QABosonRules[{{a1, c1}, {a2, c2}}]` declare their operators noncommutative. If you type `a1 ** c1`
before that, NCAlgebra treats it as an ordinary commutative product.

## Quick start

```wolfram
SetCommutative[t, m];

bosons = QABosonRules[{a, c}];        (* a**c -> c**a + 1 *)
canon  = QACanonicalRules[{x, p}];    (* x**p -> p**x + I ħ *)

QANormalOrder[a ** a ** c, bosons]    (* c**a**a + 2 a *)
QACommutator[x, p, canon]             (* I ħ *)

(* Heisenberg picture *)
QAPropagateOperator[a, ħ ω c ** a, bosons, t]
(* E^(-I t ω) a *)

osc = p ** p/(2 m) + m ω^2 x ** x/2;
QAHeisenbergSolve[{x, p}, osc, canon, t, "Assumptions" -> ħ > 0 && ω > 0 && m > 0]
(* {x Cos[t ω] + p Sin[t ω]/(m ω), p Cos[t ω] - m ω x Sin[t ω]} *)

(* Expectation values *)
QAFockExpectation[c ** c ** a ** a, {a, c}, n]     (* (n - 1) n *)
QACoherentExpectation[c ** a, {a, c}, α]           (* α Conjugate[α] *)

(* Check numerically *)
QAToMatrix[c ** a, QABosonMatrices[{a, c}, 3]] // MatrixForm
```

See [`Examples/Examples.wl`](Examples/Examples.wl) for a longer tour.

## Functions

### Normal ordering and commutators

| Function | Description |
|---|---|
| `QANormalOrder[expr, rules]` | Rewrite with `rules` until nothing changes. |
| `QACommutator[x, y, rules]` | `x**y - y**x`, normal ordered. |
| `QAAntiCommutator[x, y, rules]` | `x**y + y**x`, normal ordered. |
| `QANestedCommutator[x, y, n, rules]` | `[x, [x, … [x, y]]]`, with `x` applied `n` times. |
| `QACheckRules[rules]` | Report conflicting, cyclic or non-confluent rules. |

### Ready-made algebras

| Function | Rules |
|---|---|
| `QABosonRules[{a, c}]`, `QABosonRules[{{a1, c1}, …}]` | `[a, c] = 1`; different modes commute |
| `QAFermionRules[{f, g}]`, multi-mode | `{f, g} = 1`, `f² = g² = 0`; modes anticommute |
| `QACanonicalRules[{x, p}]`, multi-mode | `[x, p] = iħ` |
| `QASpinRules[{jx, jy, jz}]` | `[jx, jy] = iħ jz` and cyclic permutations |
| `QASpinLadderRules[{jp, jm, jz}]` | `[jz, j±] = ±ħ j±`, `[jp, jm] = 2ħ jz` |
| `QAPauliRules[{sx, sy, sz}]` | `σᵢσⱼ = δᵢⱼ + i εᵢⱼₖ σₖ` |
| `QACommutingRules[{group1, group2, …}]` | Operators in different groups commute |

Combine rule sets with `Join`. For example, for a Jaynes–Cummings model:

```wolfram
jc = Join[QABosonRules[{a, c}], QAPauliRules[{sx, sy, sz}], QACommutingRules[{{a, c}, {sx, sy, sz}}]];
```

### Adjoints

| Function | Description |
|---|---|
| `QADagger[expr, {a, c}]` | Hermitian conjugate. Scalars are treated as real unless listed in `"ComplexParameters"`. |
| `QAHermitianQ[expr, {a, c}, rules]` | Is `expr` Hermitian? |

### Dynamics and transformations

| Function | Description |
|---|---|
| `QAHeisenbergEOM[op, H, rules]` | `(i/ħ)[H, op]` |
| `QAPropagateOperator[op, H, rules, t]` | `e^{iHt/ħ} op e^{-iHt/ħ}`; optional fifth argument `t0` |
| `QAHeisenbergSolve[{ops…}, H, rules, t]` | Solve the equations of motion exactly (linear closure + `MatrixExp`) |
| `QAExpandDualExponential[op, A, rules]` | `e^A op e^-A` |
| `QABCH[x, y, rules, n]` | `Z` with `e^x e^y = e^Z`, through order `n` (default 4) |
| `QAOperatorSubstitute[expr, {a -> …}, rules]` | Substitute operators in every product, then normal order |
| `QADisplace[expr, {a, c}, α]` | `D(α)† expr D(α)` |
| `QASqueeze[expr, {a, c}, r, θ]` | `S(ξ)† expr S(ξ)`, with `ξ = r e^{iθ}` |

### Expectation values (bosonic modes)

| Function | Description |
|---|---|
| `QAFockExpectation[expr, {a, c}, n]` | `⟨n\| expr \|n⟩` |
| `QAFockMatrixElement[expr, {a, c}, m, n]` | `⟨m\| expr \|n⟩` |
| `QACoherentExpectation[expr, {a, c}, α]` | `⟨α\| expr \|α⟩` |
| `QAThermalExpectation[expr, {a, c}, nbar]` | Thermal state with mean occupation `nbar` |

All of these accept multi-mode lists, e.g. `QAFockExpectation[expr, {{a1, c1}, {a2, c2}}, {n1, n2}]`.

### Matrix representations

| Function | Description |
|---|---|
| `QAToMatrix[expr, reps]` | Matrix of `expr`. `reps` is either one list of rules, or a list of them (one per subsystem, combined as a tensor product). |
| `QABosonMatrices[{a, c}, nmax]` | Truncated Fock-space matrices (multi-mode supported). |
| `QASpinMatrices[{jx, jy, jz}, j]` | Spin-`j` matrices. |
| `QASpinLadderMatrices[{jp, jm, jz}, j]` | Spin-`j` ladder matrices. |
| `QAPauliMatrices[{sx, sy, sz}]` | Pauli matrices. |

### Ordering prescriptions, angular momentum and display

| Function | Description |
|---|---|
| `QASymmetrize[expr, a, b]` | Average of the "all `a` first" and "all `b` first" orderings. |
| `QAWeylOrder[expr, {a, b, …}]` | Average over all distinct orderings. |
| `QACGTable[j1, j2, M]`, `QACGTables[j1, j2]` | Clebsch–Gordan tables. |
| `QACoupledState[j1, j2, J, M]` | `\|J M⟩` in the product basis. |
| `QASetDaggerFormat[{a, c}]` | Display `c` as `a†` (undo with `QAClearDaggerFormat`). |

## Options

| Option | Default | Used by |
|---|---|---|
| `"MaxIterations"` | `1000` | Normal ordering |
| `"Method"` | `Automatic` | Normal ordering: `"WordRewriting"` or `"NCReplace"` |
| `"Expansion"` | `Automatic` | `e^A B e^-A`: `"Series"`, `"LinearClosure"`, or the series with linear closure as a fallback |
| `"MaxOrder"` | `20` | Series length |
| `"ClosedForm"` | `True` | Sum the series exactly when `Cₙ = μ Cₙ₋ₚ` |
| `"MaxPeriod"` | `2` | Largest period `p` checked |
| `"MaxDimension"` | `20` | Largest operator space tried by linear closure |
| `"Assumptions"` | `$Assumptions` | Used to simplify the coefficients in the result |
| `"HBar"` | `Global`ħ` | Use `1` for natural units |
| `"ComplexParameters"` | `{}` | `QADagger`: scalars to conjugate |

## How it works

**Normal ordering.** When every rule has the form *word → polynomial* (all the ready-made algebras
do), QuantumAlgebra uses its own rewriting engine. It treats `a^2 ** c` as `a ** a ** c`, so powers
are handled correctly. Rules with patterns fall back to NCAlgebra's `NCReplaceRepeated`.

**Exact Heisenberg dynamics.** `QAPropagateOperator` sums `B + λ[A,B] + λ²/2! [A,[A,B]] + …`.
It stops in one of four ways:

1. A nested commutator vanishes, so the series is finite.
2. A commutator repeats up to a scalar factor (`Cₙ = μ Cₙ₋ₚ`), so the remaining terms are summed exactly.
3. The words reached by repeated commutators span a finite space. The equations of motion are then
   linear, and the result is `MatrixExp` of that system.
4. Otherwise you get the truncated series and a `::trunc` warning.

**Consistency.** `QACheckRules` checks every overlap between rule left-hand sides
(Knuth–Bendix critical pairs). If two ways of reducing an overlap give different results,
the rules don't define a unique normal form. That usually means a missing rule or a
relation that violates the Jacobi identity.

## Notes

- `QASymmetrize` and `QAWeylOrder` are ordering *prescriptions*: they ignore commutation relations.
  They only reorder adjacent runs of the listed operators.
- `QABosonMatrices` truncates the Fock space, so `[a, c] = 1` fails in the highest level.
- Rule sets for different fermion species built separately commute with each other. Use one
  multi-mode `QAFermionRules` call if they should anticommute.

## Tests

```wolfram
TestReport["/path/to/QuantumAlgebra/Tests/QuantumAlgebra.wlt"]
```

## License

MIT. See [LICENSE](LICENSE).

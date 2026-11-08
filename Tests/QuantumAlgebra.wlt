(* Run with:  TestReport["Tests/QuantumAlgebra.wlt"]  from the repository root. *)

BeginTestSection["QuantumAlgebra"]

VerificationTest[
  Needs["QuantumAlgebra`"];
  (* single lowercase letters are noncommutative by default in NCAlgebra *)
  SetCommutative[t, m, n, r];
  (* rule sets declare their (multi-letter) operators noncommutative, so build them first *)
  bosons  = QABosonRules[{a, c}];
  bosons2 = QABosonRules[{{a1, c1}, {a2, c2}}];
  fermions = QAFermionRules[{f, g}];
  canon   = QACanonicalRules[{x, p}];
  spin    = QASpinRules[{jx, jy, jz}];
  ladder  = QASpinLadderRules[{jp, jm, jz}];
  pauli   = QAPauliRules[{sx, sy, sz}];
  osc     = p ** p/(2 m) + m \[Omega]^2 x ** x/2;
  pos     = \[HBar] > 0 && \[Omega] > 0 && m > 0 && t \[Element] Reals;
  MemberQ[$Packages, "QuantumAlgebra`"],
  True,
  TestID -> "Load"
]


(* ---- Normal ordering & commutators ---- *)

VerificationTest[
  QACommutator[a, c, a ** c -> c ** a + 1],
  1,
  TestID -> "Commutator-ladder-handwritten"
]

VerificationTest[
  NCExpand[QANormalOrder[a ** a ** c, bosons] - (c ** a ** a + 2 a)],
  0,
  TestID -> "NormalOrder-aac"
]

VerificationTest[
  NCExpand[QANormalOrder[a ** a ** c, bosons, "Method" -> "NCReplace"] - (c ** a ** a + 2 a)],
  0,
  TestID -> "NormalOrder-NCReplace-fallback"
]

VerificationTest[
  QACommutator[x, p, canon],
  I \[HBar],
  TestID -> "Commutator-xp"
]

VerificationTest[
  NCExpand[QAAntiCommutator[a, c, bosons] - (2 c ** a + 1)],
  0,
  TestID -> "AntiCommutator"
]

VerificationTest[
  QANestedCommutator[x ** x, p, 2, canon],
  0,
  TestID -> "NestedCommutator"
]


(* ---- Ready-made algebras ---- *)

VerificationTest[
  {QACommutator[a1, c1, bosons2], QACommutator[a1, c2, bosons2], NCExpand[QANormalOrder[a2 ** a1, bosons2] - a1 ** a2]},
  {1, 0, 0},
  TestID -> "Bosons-multimode"
]

VerificationTest[
  {QAAntiCommutator[f, g, fermions], QANormalOrder[g ** g, fermions], QANormalOrder[f ** g ** g, fermions]},
  {1, 0, 0},
  TestID -> "Fermions"
]

VerificationTest[
  {QACommutator[jx, jy, spin], QACommutator[jp, jm, ladder], QACommutator[jz, jp, ladder]},
  {I \[HBar] jz, 2 \[HBar] jz, \[HBar] jp},
  SameTest -> ((NCExpand /@ (#1 - #2)) === {0, 0, 0} &),
  TestID -> "Spin"
]

VerificationTest[
  QANormalOrder[sx ** sy ** sz, pauli],
  I,
  TestID -> "Pauli"
]

VerificationTest[
  Lookup[QACheckRules[#], "Confluent"] & /@ {bosons, bosons2, fermions, spin, ladder, pauli},
  {True, True, True, True, True, True},
  TestID -> "CheckRules-presets-confluent"
]

VerificationTest[
  QACheckRules[{b ** a -> 1, c ** b -> 1}]["Confluent"],
  False,
  TestID -> "CheckRules-nonconfluent"
]


(* ---- Adjoints ---- *)

VerificationTest[
  NCExpand[QADagger[2 I a ** c ** c, {a, c}] - (-2 I a ** a ** c)],
  0,
  TestID -> "Dagger"
]

VerificationTest[
  NCExpand[QADagger[\[Alpha] c, {a, c}, "ComplexParameters" -> {\[Alpha]}] - Conjugate[\[Alpha]] a],
  0,
  TestID -> "Dagger-complex-parameter"
]

VerificationTest[
  {QAHermitianQ[c ** a, {a, c}, bosons], QAHermitianQ[I c ** a, {a, c}, bosons]},
  {True, False},
  TestID -> "HermitianQ"
]


(* ---- Similarity transforms ---- *)

VerificationTest[
  NCExpand[QAExpandDualExponential[x, \[Kappa] p, canon] - (x - I \[Kappa] \[HBar])],
  0,
  TestID -> "DualExp-terminating"
]

VerificationTest[
  NCExpand[QAExpandDualExponential[a, \[Kappa] c ** a, bosons, "Expansion" -> "LinearClosure"] - Exp[-\[Kappa]] a],
  0,
  TestID -> "DualExp-linear-closure"
]

VerificationTest[
  QAExpandDualExponential[b, a, {}, "MaxOrder" -> 3, "MaxDimension" -> 5];,
  Null,
  {QAExpandDualExponential::trunc},
  TestID -> "DualExp-truncation-warning"
]


(* ---- Dynamics ---- *)

VerificationTest[
  QAHeisenbergEOM[a, \[HBar] \[Omega] c ** a, bosons],
  -I \[Omega] a,
  SameTest -> (NCExpand[#1 - #2] === 0 &),
  TestID -> "HeisenbergEOM"
]

VerificationTest[
  Simplify[QAPropagateOperator[a, \[HBar] \[Omega] c ** a, bosons, t] / a],
  Exp[-I \[Omega] t],
  TestID -> "Propagate-ladder"
]

VerificationTest[
  Module[{res = QAPropagateOperator[x, osc, canon, t, "Assumptions" -> pos]},
    Simplify[{Coefficient[res, x], Coefficient[res, p]}, pos]],
  {Cos[\[Omega] t], Sin[\[Omega] t]/(m \[Omega])},
  TestID -> "Propagate-oscillator-series"
]

VerificationTest[
  Module[{res = QAHeisenbergSolve[{x, p}, osc, canon, t, "Assumptions" -> pos]},
    Simplify[{Coefficient[res[[1]], x], Coefficient[res[[1]], p], Coefficient[res[[2]], x]}, pos]],
  {Cos[\[Omega] t], Sin[\[Omega] t]/(m \[Omega]), -m \[Omega] Sin[\[Omega] t]},
  TestID -> "HeisenbergSolve-oscillator"
]

VerificationTest[
  Simplify[QAPropagateOperator[a, \[HBar] \[Omega] c ** a, bosons, t, t0] / a],
  Exp[-I \[Omega] (t - t0)],
  SameTest -> (Simplify[#1 - #2] === 0 &),
  TestID -> "Propagate-t0"
]


(* ---- BCH and substitutions ---- *)

VerificationTest[
  NCExpand[QABCH[\[Alpha] c, -Conjugate[\[Alpha]] a, bosons] -
    (\[Alpha] c - Conjugate[\[Alpha]] a + \[Alpha] Conjugate[\[Alpha]]/2)],
  0,
  TestID -> "BCH-displacement"
]

VerificationTest[
  NCExpand[QADisplace[c ** a, {a, c}, \[Alpha]] -
    (c ** a + \[Alpha] c + Conjugate[\[Alpha]] a + \[Alpha] Conjugate[\[Alpha]])],
  0,
  TestID -> "Displace"
]

VerificationTest[
  NCExpand[QASqueeze[a, {a, c}, r] - (a Cosh[r] - c Sinh[r])],
  0,
  TestID -> "Squeeze"
]


(* ---- Expectation values ---- *)

VerificationTest[
  {QAFockExpectation[c ** a, {a, c}, n], QAFockExpectation[a ** c, {a, c}, n],
   QAFockExpectation[c ** c ** a ** a, {a, c}, n]},
  {n, n + 1, n (n - 1)},
  SameTest -> (Simplify[#1 - #2] === {0, 0, 0} &),
  TestID -> "FockExpectation"
]

VerificationTest[
  QAFockMatrixElement[c, {a, c}, n + 1, n],
  Sqrt[n + 1],
  TestID -> "FockMatrixElement"
]

VerificationTest[
  Simplify[QACoherentExpectation[c ** a, {a, c}, \[Alpha]] - \[Alpha] Conjugate[\[Alpha]]],
  0,
  TestID -> "CoherentExpectation"
]

VerificationTest[
  QAThermalExpectation[a ** c, {a, c}, nb],
  nb + 1,
  TestID -> "ThermalExpectation"
]


(* ---- Matrices ---- *)

VerificationTest[
  Normal @ QAToMatrix[c ** a, QABosonMatrices[{a, c}, 3]],
  DiagonalMatrix[{0, 1, 2, 3}],
  TestID -> "Matrix-number-operator"
]

VerificationTest[
  Simplify @ Normal @ QAToMatrix[jx ** jy - jy ** jx - I \[HBar] jz, QASpinMatrices[{jx, jy, jz}, 1]],
  ConstantArray[0, {3, 3}],
  TestID -> "Matrix-spin-1-commutator"
]

VerificationTest[
  Normal @ QAToMatrix[c ** a ** sz, {QABosonMatrices[{a, c}, 1], QAPauliMatrices[{sx, sy, sz}]}],
  KroneckerProduct[{{0, 0}, {0, 1}}, PauliMatrix[3]],
  TestID -> "Matrix-tensor-product"
]


(* ---- Ordering prescriptions ---- *)

VerificationTest[
  NCExpand[QASymmetrize[x ** p, x, p] - (x ** p + p ** x)/2],
  0,
  TestID -> "Symmetrize-xp"
]

VerificationTest[
  NCExpand[QAWeylOrder[x ** x ** p, {x, p}] - (x ** x ** p + x ** p ** x + p ** x ** x)/3],
  0,
  TestID -> "Weyl-xxp"
]


(* ---- Angular momentum coupling ---- *)

VerificationTest[
  First[QACGTable[1/2, 1/2, 0]],
  {{1/Sqrt[2], 1/Sqrt[2]}, {1/Sqrt[2], -1/Sqrt[2]}},
  TestID -> "CG-two-spin-half"
]

VerificationTest[
  QACGTable[1/3, 1, 0],
  $Failed,
  {QACGTable::args},
  TestID -> "CG-bad-args"
]

VerificationTest[
  QACoupledState[1/2, 1/2, 0, 0],
  <|{1/2, -1/2} -> 1/Sqrt[2], {-1/2, 1/2} -> -1/Sqrt[2]|>,
  TestID -> "CoupledState-singlet"
]

EndTestSection[]

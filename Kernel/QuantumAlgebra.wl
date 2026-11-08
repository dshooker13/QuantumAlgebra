(* ::Package:: *)

(* :Title:    QuantumAlgebra                                                     *)
(* :Author:   Daniel Hooker <dshooker13@gmail.com>                               *)
(* :Summary:  Operator algebra for quantum mechanics, built on NCAlgebra:         *)
(*            normal ordering, ready-made algebras, commutators, Heisenberg       *)
(*            picture, BCH, expectation values, matrix representations,           *)
(*            ordering prescriptions and Clebsch-Gordan tools.                    *)
(* :Package Version: 1.0.1                                                      *)
(* :Mathematica Version: 12.1+                                                   *)
(* :Requirements: NCAlgebra 5.0+ (https://github.com/NCAlgebra/NC)               *)
(* :License: MIT                                                                 *)


(* ::Section:: *)
(*Load NCAlgebra*)


(*
   NCAlgebra is loaded *before* BeginPackage so that its sub-contexts
   (NonCommutativeMultiply`, NCReplace`, ...) end up on the user's $ContextPath.
   Older NCAlgebra installs need NC` loaded first to set up paths; newer ones
   don't have it, so failure there is silenced.
*)
If[!MemberQ[$Packages, "NCAlgebra`"],
  Quiet[Needs["NC`"]];
  Needs["NCAlgebra`"]
];

QuantumAlgebra`Private`$ncContexts =
  Select[$ContextPath, StringStartsQ[#, "NC"] || # === "NonCommutativeMultiply`" &];


BeginPackage["QuantumAlgebra`"];


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Normal ordering and commutators*)


QANormalOrder::usage =
"QANormalOrder[expr, rules] rewrites expr with rules until nothing changes, putting it in the order the rules define.
Example: QANormalOrder[a**a**c, a**c -> c**a + 1] gives c**a**a + 2 a.
When every rule has the form word -> polynomial (e.g. a**c -> c**a + 1, f**f -> 0), a built-in word-rewriting \
engine is used; it handles powers such as a^2**c correctly. Rules with patterns fall back to NCReplaceRepeated.
Options: \"MaxIterations\" (1000), \"Method\" (Automatic, \"WordRewriting\" or \"NCReplace\").";

QACommutator::usage =
"QACommutator[x, y, rules] gives the commutator x**y - y**x, normal ordered with rules.";

QAAntiCommutator::usage =
"QAAntiCommutator[x, y, rules] gives the anticommutator x**y + y**x, normal ordered with rules.";

QANestedCommutator::usage =
"QANestedCommutator[x, y, n, rules] gives [x, [x, ... [x, y]]] with x applied n times, normal ordered with rules.";

QACheckRules::usage =
"QACheckRules[rules] checks a set of word -> polynomial rules for conflicting rules, rewrite cycles and \
non-confluence (ambiguous overlaps whose two reductions differ, e.g. relations that violate the Jacobi identity). \
Returns an Association whose \"Confluent\" key is True when the rules define a consistent normal form.";


(* ::Subsection:: *)
(*Ready-made algebras*)


QABosonRules::usage =
"QABosonRules[{a, c}] gives the normal-ordering rule a**c -> c**a + 1 for a bosonic mode with annihilator a and creator c.
QABosonRules[{{a1, c1}, {a2, c2}, ...}] gives rules for several modes: [ai, cj] = \[Delta]ij and operators of \
different modes commute. Creators end up to the left of annihilators.
The operators are declared noncommutative automatically.";

QAFermionRules::usage =
"QAFermionRules[{f, g}] gives normal-ordering rules for a fermionic mode with annihilator f and creator g: \
{f, g} = 1 and f**f = g**g = 0.
QAFermionRules[{{f1, g1}, {f2, g2}, ...}] gives rules for several modes, which mutually anticommute.";

QACanonicalRules::usage =
"QACanonicalRules[{x, p}] gives x**p -> p**x + I \[HBar].
QACanonicalRules[{{x1, p1}, {x2, p2}, ...}] gives rules for several degrees of freedom; different ones commute.
Option \"HBar\" (default Global`\[HBar]).";

QASpinRules::usage =
"QASpinRules[{jx, jy, jz}] gives angular-momentum rules [ji, jj] = I \[HBar] \[Epsilon]ijk jk, ordering \
products as jx, jy, jz. Option \"HBar\" (default Global`\[HBar]).";

QASpinLadderRules::usage =
"QASpinLadderRules[{jp, jm, jz}] gives angular-momentum rules in the ladder basis, [jz, j\[PlusMinus]] = \
\[PlusMinus]\[HBar] j\[PlusMinus] and [jp, jm] = 2 \[HBar] jz, ordering products as jp, jz, jm. \
Option \"HBar\" (default Global`\[HBar]).";

QAPauliRules::usage =
"QAPauliRules[{sx, sy, sz}] gives the Pauli-matrix product rules \[Sigma]i \[Sigma]j = \[Delta]ij + I \
\[Epsilon]ijk \[Sigma]k.";

QACommutingRules::usage =
"QACommutingRules[{group1, group2, ...}] gives rules making operators from different groups commute, moving \
earlier groups to the left. Use it to combine algebras, e.g. Join[QABosonRules[{a, c}], QAPauliRules[{sx, sy, sz}], \
QACommutingRules[{{a, c}, {sx, sy, sz}}]].";


(* ::Subsection:: *)
(*Adjoints*)


QADagger::usage =
"QADagger[expr, {a, c}] gives the Hermitian conjugate of expr: products are reversed, a and c are swapped, \
I -> -I, and every other operator is treated as self-adjoint.
QADagger[expr, {{a1, c1}, {a2, c2}, ...}] swaps several pairs. QADagger[expr] treats every operator as self-adjoint.
Scalar symbols are treated as real unless listed in the option \"ComplexParameters\".";

QAHermitianQ::usage =
"QAHermitianQ[expr, pairs, rules] gives True if expr equals its Hermitian conjugate after normal ordering with rules.";


(* ::Subsection:: *)
(*Transformations and dynamics*)


QAExpandDualExponential::usage =
"QAExpandDualExponential[op, a, rules] gives Exp[a] ** op ** Exp[-a] = op + [a, op] + 1/2! [a, [a, op]] + ..., \
normal ordering each commutator with rules.
The result is exact when the series terminates, when a commutator repeats up to a scalar factor, or when the \
operators generated close on a finite-dimensional linear space (solved with MatrixExp).
Options: \"Expansion\" (Automatic, \"Series\" or \"LinearClosure\"), \"MaxOrder\" (20), \"ClosedForm\" (True), \
\"MaxPeriod\" (2), \"MaxDimension\" (20), \"Assumptions\" ($Assumptions), plus the QANormalOrder options.";

QAPropagateOperator::usage =
"QAPropagateOperator[op, H, rules, t] gives the Heisenberg-picture operator \
Exp[I H t/\[HBar]] ** op ** Exp[-I H t/\[HBar]] for a time-independent Hamiltonian H.
QAPropagateOperator[op, H, rules, t, t0] uses the elapsed time t - t0.
Takes the options of QAExpandDualExponential plus \"HBar\" (default Global`\[HBar]).";

QAPropogateOperator::usage =
"QAPropogateOperator is a deprecated misspelling of QAPropagateOperator, kept for backward compatibility.";

QAHeisenbergEOM::usage =
"QAHeisenbergEOM[op, H, rules] gives the Heisenberg equation of motion d op/dt = (I/\[HBar]) [H, op], \
normal ordered with rules. Option \"HBar\" (default Global`\[HBar]).";

QAHeisenbergSolve::usage =
"QAHeisenbergSolve[{op1, op2, ...}, H, rules, t] solves the Heisenberg equations of motion exactly when the \
operators generated by repeated commutators with H span a finite-dimensional space, and returns {op1(t), op2(t), ...}.
QAHeisenbergSolve[op, H, rules, t] returns a single operator. Returns $Failed if no closure is found within \
\"MaxDimension\" (default 20) words.";

QABCH::usage =
"QABCH[x, y, rules, n] gives Z with Exp[x] ** Exp[y] = Exp[Z] (Baker-Campbell-Hausdorff), through order n in x \
and y (default 4), normal ordered with rules. The result is exact when the higher commutators vanish, e.g. when \
[x, y] is a scalar.";

QAOperatorSubstitute::usage =
"QAOperatorSubstitute[expr, {a -> e1, b -> e2, ...}, rules] replaces each operator in every product of expr by \
the given expression, then normal orders with rules. Use it for any linear (Bogoliubov-type) transformation.";

QADisplace::usage =
"QADisplace[expr, {a, c}, \[Alpha]] gives D(\[Alpha])\[ConjugateTranspose] ** expr ** D(\[Alpha]) for the \
displacement operator D(\[Alpha]), i.e. a -> a + \[Alpha], c -> c + Conjugate[\[Alpha]], normal ordered.
QADisplace[expr, {{a1, c1}, ...}, {\[Alpha]1, ...}] displaces several modes.";

QASqueeze::usage =
"QASqueeze[expr, {a, c}, r, \[Theta]] gives S(\[Xi])\[ConjugateTranspose] ** expr ** S(\[Xi]) for the \
squeeze operator S(\[Xi]) = Exp[(Conjugate[\[Xi]] a^2 - \[Xi] c^2)/2], \[Xi] = r Exp[I \[Theta]]: \
a -> a Cosh[r] - c Exp[I \[Theta]] Sinh[r]. \[Theta] defaults to 0.";


(* ::Subsection:: *)
(*Expectation values*)


QAFockExpectation::usage =
"QAFockExpectation[expr, {a, c}, n] gives <n| expr |n> for a bosonic Fock state.
QAFockExpectation[expr, {{a1, c1}, ...}, {n1, ...}] uses a multi-mode Fock state.";

QAFockMatrixElement::usage =
"QAFockMatrixElement[expr, {a, c}, m, n] gives <m| expr |n> between bosonic Fock states. Multi-mode forms \
take lists of modes and occupation numbers.";

QACoherentExpectation::usage =
"QACoherentExpectation[expr, {a, c}, \[Alpha]] gives <\[Alpha]| expr |\[Alpha]> for a coherent state. \
Multi-mode forms take lists of modes and amplitudes.";

QAThermalExpectation::usage =
"QAThermalExpectation[expr, {a, c}, nbar] gives the thermal-state expectation value with mean occupation nbar. \
Multi-mode forms take lists of modes and occupations.";


(* ::Subsection:: *)
(*Matrix representations*)


QAToMatrix::usage =
"QAToMatrix[expr, {op1 -> m1, ...}] gives the matrix of expr, with each operator replaced by its matrix.
QAToMatrix[expr, {reps1, reps2, ...}] treats each list of rules as a separate subsystem and embeds the operators \
in the tensor product, in the order given.";

QABosonMatrices::usage =
"QABosonMatrices[{a, c}, nmax] gives {a -> matrix, c -> matrix} on the Fock space truncated at nmax quanta.
QABosonMatrices[{{a1, c1}, ...}, nmax] gives matrices on the tensor product of all modes (nmax may be a list). \
Note that truncation breaks [a, c] = 1 in the highest level.";

QASpinMatrices::usage =
"QASpinMatrices[{jx, jy, jz}, j] gives the spin-j matrices, in the basis m = j, j-1, ..., -j. Option \"HBar\".";

QASpinLadderMatrices::usage =
"QASpinLadderMatrices[{jp, jm, jz}, j] gives the spin-j ladder and jz matrices. Option \"HBar\".";

QAPauliMatrices::usage =
"QAPauliMatrices[{sx, sy, sz}] gives the Pauli matrices.";


(* ::Subsection:: *)
(*Ordering prescriptions*)


QASymmetrize::usage =
"QASymmetrize[expr, a, b] gives (expr_ab + expr_ba)/2, where expr_ab puts every a to the left of every \
b in each product and expr_ba does the reverse. Commutation relations are ignored, so \
QASymmetrize[x**x**p, x, p] gives (x**x**p + p**x**x)/2. Only adjacent runs of a and b are \
reordered; other operators stay in place.";

QAWeylOrder::usage =
"QAWeylOrder[expr, {a, b, ...}] Weyl-orders (fully symmetrizes) every product in expr with respect \
to the listed operators: each word is replaced by the average of all its distinct orderings. Commutation \
relations are ignored, so QAWeylOrder[x**x**p, {x, p}] gives (x**x**p + x**p**x + p**x**x)/3. \
Only adjacent runs of the listed operators are reordered; other operators stay in place.";


(* ::Subsection:: *)
(*Angular momentum coupling*)


QACGTable::usage =
"QACGTable[j1, j2, M] gives the table of Clebsch-Gordan coefficients <j1 m1; j2 m2 | J M> for fixed M. \
Rows are labeled by (m1, m2), columns by J.";

QACGTables::usage =
"QACGTables[j1, j2] gives the Clebsch-Gordan tables for every M from j1 + j2 down to -(j1 + j2).";

QACoupledState::usage =
"QACoupledState[j1, j2, J, M] gives the coupled state |J M> in the product basis, as an Association \
{m1, m2} -> <j1 m1; j2 m2 | J M> of nonzero coefficients.";


(* ::Subsection:: *)
(*Display*)


QASetDaggerFormat::usage =
"QASetDaggerFormat[{a, c}] makes c display as a\[Dagger]. QASetDaggerFormat[{{a1, c1}, ...}] does this for several pairs.";

QAClearDaggerFormat::usage =
"QAClearDaggerFormat[{a, c}] undoes QASetDaggerFormat.";


(* ::Section:: *)
(*Messages*)


General::qancparam =
"`1` is noncommutative, but is used here as a scalar. NCAlgebra treats every single lowercase letter as \
noncommutative by default; declare it with SetCommutative[`1`].";
General::qamodes =
"Modes must be given as {a, c} or {{a1, c1}, {a2, c2}, ...}.";
General::qanonpoly =
"The expression must be a polynomial in the operators (sums of scalar multiples of products).";
General::qalength =
"Expected one value per mode.";

QANormalOrder::maxit =
"Normal ordering did not converge within `1` iterations; returning the last result. \
Check the rules for cycles (QACheckRules can help) or raise \"MaxIterations\".";
QANormalOrder::norw =
"The rules are not all of the form word -> polynomial; using NCReplaceRepeated instead.";

QACheckRules::simple =
"QACheckRules only handles rules of the form word -> polynomial, without patterns.";

QAExpandDualExponential::trunc =
"The nested-commutator series did not close; truncated after order `1`. \
Raise \"MaxOrder\" or \"MaxDimension\", or check that the series closes.";
QAPropagateOperator::trunc = QAExpandDualExponential::trunc;

QAExpandDualExponential::closure =
"The operators generated by repeated commutators did not close on a space of at most `1` words.";
QAPropagateOperator::closure = QAExpandDualExponential::closure;
QAHeisenbergSolve::closure = QAExpandDualExponential::closure;

QAPropogateOperator::deprec =
"QAPropogateOperator is deprecated; use QAPropagateOperator.";

QAFockExpectation::ops = QAFockMatrixElement::ops = QACoherentExpectation::ops = QAThermalExpectation::ops =
"The expression contains operators other than the given bosonic modes.";

QAToMatrix::reps = "Representations must be a list of rules or a list of lists of rules.";
QAToMatrix::ops = "No matrix given for the operators `1`.";

QACGTable::args =
"j1 and j2 must be non-negative integers or half-integers, and M must be one of j1 + j2, j1 + j2 - 1, ..., -(j1 + j2).";
QACGTables::args =
"j1 and j2 must be non-negative integers or half-integers.";
QACoupledState::args =
"Need non-negative (half-)integers j1, j2 with |j1 - j2| <= J <= j1 + j2 and M in -J, ..., J.";


Begin["`Private`"];

$ContextPath = DeleteDuplicates[Join[$ContextPath, $ncContexts]];


(* ::Section:: *)
(*Internal helpers*)


(* ::Subsection:: *)
(*Scalars, words and coefficient tables*)


(* Zero test that survives scalar coefficients that only cancel after simplification. *)
iZeroQ[x_] := Module[{parts},
  If[x === 0, Return[True]];
  parts = iSplitTerm /@ iTerms[NCExpand[x]];
  If[MemberQ[parts, $Failed], Return[False]];
  AllTrue[
    Values @ Merge[(#[[2]] -> #[[1]]) & /@ parts, Total],
    TrueQ[Quiet[Simplify[Expand[#]]] === 0] &
  ]
];

(* Cheap zero test used inside loops. *)
iNonzeroCoefQ[k_] := !(k === 0 || TrueQ[Expand[k] === 0]);

(* Is s a noncommutative atom (symbol or subscripted symbol)? *)
iNCAtomQ[s_Symbol] := !TrueQ[CommutativeQ[s]];
iNCAtomQ[s : Subscript[_, __]] := !TrueQ[CommutativeQ[s]];
iNCAtomQ[_] := False;

(* Does x contain no noncommutative operators at all? *)
iScalarQ[x_] := FreeQ[x, _NonCommutativeMultiply] && FreeQ[x, _?iNCAtomQ];

(* Warn when a parameter meant as a scalar is noncommutative (e.g. a bare t). *)
iCheckScalars[caller_Symbol, xs_List] :=
  Scan[If[!iScalarQ[#], Message[MessageName[caller, "qancparam"], #]] &, xs];

iTerms[e_Plus] := List @@ e;
iTerms[e_] := {e};

(* Flatten a product into its list of letters: a**b^2**c -> {a, b, b, c}. *)
iLetters[w_NonCommutativeMultiply] := Join @@ (iLetters /@ List @@ w);
iLetters[Power[x_, n_Integer?Positive]] /; !iScalarQ[x] := Join @@ ConstantArray[iLetters[x], n];
iLetters[x_] := {x};

iWord[{}] := 1;
iWord[{x_}] := x;
iWord[l_List] := NonCommutativeMultiply @@ l;

(* Split one expanded term into {scalar coefficient, letters}; $Failed if it isn't a plain monomial. *)
iSplitTerm[t_] := Module[{factors, nc},
  factors = If[Head[t] === Times, List @@ t, {t}];
  nc = Select[factors, !iScalarQ[#] &];
  Which[
    nc === {}, {t, {}},
    Length[nc] == 1, {Times @@ Select[factors, iScalarQ], iLetters[First[nc]]},
    True, $Failed
  ]
];

(* <|letters -> coefficient|> for an expression, dropping zero coefficients; $Failed if not a polynomial. *)
iWordCoefficients[e_] := Module[{parts},
  parts = iSplitTerm /@ iTerms[NCExpand[e]];
  If[MemberQ[parts, $Failed],
    $Failed,
    Select[Merge[(#[[2]] -> #[[1]]) & /@ parts, Total], !iZeroQ[#] &]
  ]
];

iFromCoefficients[assoc_Association] := Total[KeyValueMap[#2 iWord[#1] &, assoc]];

(* If x == mu y for a nonzero scalar mu, return mu; otherwise None. *)
iProportionality[x_, y_] := Module[{cx = iWordCoefficients[x], cy = iWordCoefficients[y], keys, mu},
  If[cx === $Failed || cy === $Failed || cy === <||> || Sort[Keys[cx]] =!= Sort[Keys[cy]],
    None,
    keys = Keys[cy];
    mu = Simplify[cx[First[keys]] / cy[First[keys]]];
    If[iScalarQ[mu] && AllTrue[keys, iZeroQ[cx[#] - mu cy[#]] &], mu, None]
  ]
];

(* Pick the simpler of the exponential and trigonometric forms of a coefficient. *)
iSimplifyCoef[x_, assum_] := Module[{s1, s2},
  s1 = Quiet @ Simplify[x, Assumptions -> assum];
  s2 = Quiet @ Simplify[ExpToTrig[x], Assumptions -> assum];
  If[LeafCount[s2] < LeafCount[s1], s2, s1]
];

iList[x_List] := x;
iList[x_] := {x};


(* ::Subsection:: *)
(*Modes*)


iModes[m : {Except[_List], Except[_List]}] := {m};
iModes[m : {{_, _} ..}] := m;
iModes[_] := $Failed;

iSetNC[ops_List] := Quiet[SetNonCommutative @@ ops];


(* ::Subsection:: *)
(*Word-rewriting engine*)


iPatternFreeQ[e_] := FreeQ[e,
  _Pattern | _Blank | _BlankSequence | _BlankNullSequence | _Condition | _PatternTest |
  _Alternatives | _Optional | _Repeated | _RepeatedNull, Heads -> True];

(* lhs word -> polynomial  ==>  {lhsLetters, <|word -> coef|>}, or $Failed *)
iCompileRule[(Rule | RuleDelayed)[lhs_, rhs_]] := Module[{sl, rc},
  sl = If[iPatternFreeQ[lhs], iSplitTerm[lhs], $Failed];
  If[MatchQ[sl, {1, l_List /; Length[l] >= 2 && AllTrue[l, iNCAtomQ]}] && iPatternFreeQ[rhs],
    rc = iWordCoefficients[rhs];
    If[rc === $Failed, $Failed, {sl[[2]], rc}],
    $Failed
  ]
];
iCompileRule[_] := $Failed;

iCompileRules[rules_] := Module[{compiled = iCompileRule /@ Flatten[{rules}]},
  If[MemberQ[compiled, $Failed], $Failed, compiled]
];

(* Leftmost occurrence of any rule's left-hand side in w: {position, length, rhs} or None. *)
iFirstMatch[w_List, crules_List] := Module[{best = None, bp = Infinity, pos},
  Do[
    If[Length[rule[[1]]] <= Length[w],
      pos = SequencePosition[w, rule[[1]], 1];
      If[pos =!= {} && pos[[1, 1]] < bp,
        bp = pos[[1, 1]];
        best = {bp, Length[rule[[1]]], rule[[2]]}
      ]
    ],
    {rule, crules}
  ];
  best
];

(*
   Rewrite a coefficient table to normal form. Each pass rewrites the leftmost redex of every
   reducible word. Returns {table, convergedQ}.
*)
iRewrite[coefs_Association, crules_List, maxIter_] := Module[{cur = coefs, done = {}, next, iter = 0, m},
  While[Length[cur] > 0 && iter < maxIter,
    iter++;
    next = {};
    KeyValueMap[
      Function[{w, k},
        m = iFirstMatch[w, crules];
        If[m === None,
          AppendTo[done, w -> k],
          next = Join[next, KeyValueMap[
            Function[{rw, rk}, Join[Take[w, m[[1]] - 1], rw, Drop[w, m[[1]] + m[[2]] - 1]] -> k rk],
            m[[3]]
          ]]
        ]
      ],
      cur
    ];
    cur = Select[Merge[next, Total], iNonzeroCoefQ]
  ];
  {Select[Merge[Join[done, Normal[cur]], Total], !iZeroQ[#] &], Length[cur] == 0}
];


(* ::Section:: *)
(*Normal ordering and commutators*)


Options[QANormalOrder] = {"MaxIterations" -> 1000, "Method" -> Automatic};

QANormalOrder[expr_, rules_, OptionsPattern[]] := Module[
  {max = OptionValue["MaxIterations"], method = OptionValue["Method"], crules, wc, res},
  crules = If[method === "NCReplace", $Failed, iCompileRules[rules]];
  If[crules === $Failed && method === "WordRewriting", Message[QANormalOrder::norw]];
  wc = If[crules === $Failed, $Failed, iWordCoefficients[expr]];
  If[wc =!= $Failed,
    res = iRewrite[wc, crules, max];
    If[!res[[2]], Message[QANormalOrder::maxit, max]];
    NCExpand[iFromCoefficients[res[[1]]]],
    iNCReplaceOrder[expr, rules, max]
  ]
];

iNCReplaceOrder[expr_, rules_, max_] := Module[{steps},
  steps = FixedPointList[NCReplaceRepeated[NCExpand[#], rules] &, expr, max];
  If[Length[steps] > max && steps[[-1]] =!= steps[[-2]], Message[QANormalOrder::maxit, max]];
  Last[steps]
];

QACommutator[x_, y_, rules_, opts : OptionsPattern[QANormalOrder]] :=
  QANormalOrder[x ** y - y ** x, rules, opts];

QAAntiCommutator[x_, y_, rules_, opts : OptionsPattern[QANormalOrder]] :=
  QANormalOrder[x ** y + y ** x, rules, opts];

QANestedCommutator[x_, y_, n_Integer?NonNegative, rules_, opts : OptionsPattern[QANormalOrder]] :=
  Nest[QACommutator[x, #, rules, opts] &, y, n];


(* ::Subsection:: *)
(*Rule checking (critical pairs)*)


(* Overlaps of rule 1 with rule 2: {word, reduction via rule 1, reduction via rule 2}. *)
iCriticalPairs[{l1_, r1_}, {l2_, r2_}] := Module[{out = {}, n1 = Length[l1], n2 = Length[l2]},
  Do[
    If[Take[l1, -k] === Take[l2, k],
      AppendTo[out, {
        Join[l1, Drop[l2, k]],
        KeyMap[Join[#, Drop[l2, k]] &, r1],
        KeyMap[Join[Drop[l1, -k], #] &, r2]
      }]
    ],
    {k, 1, Min[n1, n2] - 1}
  ];
  If[n2 < n1,
    Do[
      If[l1[[s ;; s + n2 - 1]] === l2,
        AppendTo[out, {l1, r1, KeyMap[Join[Take[l1, s - 1], #, Drop[l1, s + n2 - 1]] &, r2]}]
      ],
      {s, 1, n1 - n2 + 1}
    ]
  ];
  out
];

Options[QACheckRules] = {"MaxIterations" -> 1000};

QACheckRules[rules_, OptionsPattern[]] := Module[
  {comp = iCompileRules[rules], max = OptionValue["MaxIterations"], conflicts, critical,
   nonconf = {}, nonterm = {}, p1, p2, diff},
  If[comp === $Failed,
    Message[QACheckRules::simple]; $Failed,

    conflicts = Select[GatherBy[comp, First], Length[DeleteDuplicates[#[[All, 2]]]] > 1 &][[All, 1, 1]];
    critical = Flatten[Table[iCriticalPairs[comp[[i]], comp[[k]]], {i, Length[comp]}, {k, Length[comp]}], 2];
    Do[
      p1 = iRewrite[cp[[2]], comp, max];
      p2 = iRewrite[cp[[3]], comp, max];
      If[!(p1[[2]] && p2[[2]]),
        AppendTo[nonterm, iWord[cp[[1]]]],
        diff = Select[Merge[{p1[[1]], Minus /@ p2[[1]]}, Total], !iZeroQ[#] &];
        If[diff =!= <||>, AppendTo[nonconf, iWord[cp[[1]]] -> NCExpand[iFromCoefficients[diff]]]]
      ],
      {cp, critical}
    ];
    <|
      "Confluent" -> (conflicts === {} && nonconf === {} && nonterm === {}),
      "Conflicts" -> (iWord /@ conflicts),
      "NonConfluent" -> DeleteDuplicates[nonconf],
      "NonTerminating" -> DeleteDuplicates[nonterm],
      "CriticalPairs" -> Length[critical]
    |>
  ]
];


(* ::Section:: *)
(*Ready-made algebras*)


(*
   pairs = {{r1, l1}, ...}: every l moves to the left of every r.
   r_i ** l_j -> s l_j ** r_i + delta_ij val;  same-type operators sorted by mode index (sign s).
*)
iPairRules[pairs_, s_, val_] := Module[{n = Length[pairs], r = pairs[[All, 1]], l = pairs[[All, 2]]},
  Join[
    Flatten @ Table[r[[i]] ** l[[j]] -> s l[[j]] ** r[[i]] + If[i == j, val, 0], {i, n}, {j, n}],
    Flatten @ Table[r[[j]] ** r[[i]] -> s r[[i]] ** r[[j]], {i, n}, {j, i + 1, n}],
    Flatten @ Table[l[[j]] ** l[[i]] -> s l[[i]] ** l[[j]], {i, n}, {j, i + 1, n}]
  ]
];

QABosonRules[modes_] := With[{m = iModes[modes]},
  If[m === $Failed,
    Message[QABosonRules::qamodes]; $Failed,
    iSetNC[Flatten[m]];
    iPairRules[m, 1, 1]
  ]
];

QAFermionRules[modes_] := With[{m = iModes[modes]},
  If[m === $Failed,
    Message[QAFermionRules::qamodes]; $Failed,
    iSetNC[Flatten[m]];
    Join[
      iPairRules[m, -1, 1],
      Flatten @ Table[{m[[i, 1]] ** m[[i, 1]] -> 0, m[[i, 2]] ** m[[i, 2]] -> 0}, {i, Length[m]}]
    ]
  ]
];

Options[QACanonicalRules] = {"HBar" -> Global`\[HBar]};
QACanonicalRules[modes_, OptionsPattern[]] := With[{m = iModes[modes]},
  If[m === $Failed,
    Message[QACanonicalRules::qamodes]; $Failed,
    iSetNC[Flatten[m]];
    iPairRules[m, 1, I OptionValue["HBar"]]
  ]
];

Options[QASpinRules] = {"HBar" -> Global`\[HBar]};
QASpinRules[{jx_, jy_, jz_}, OptionsPattern[]] := With[{h = OptionValue["HBar"]},
  iSetNC[{jx, jy, jz}];
  {
    jy ** jx -> jx ** jy - I h jz,
    jz ** jx -> jx ** jz + I h jy,
    jz ** jy -> jy ** jz - I h jx
  }
];

Options[QASpinLadderRules] = {"HBar" -> Global`\[HBar]};
QASpinLadderRules[{jp_, jm_, jz_}, OptionsPattern[]] := With[{h = OptionValue["HBar"]},
  iSetNC[{jp, jm, jz}];
  {
    jz ** jp -> jp ** jz + h jp,
    jm ** jp -> jp ** jm - 2 h jz,
    jm ** jz -> jz ** jm + h jm
  }
];

QAPauliRules[{sx_, sy_, sz_}] := (
  iSetNC[{sx, sy, sz}];
  {
    sx ** sx -> 1, sy ** sy -> 1, sz ** sz -> 1,
    sx ** sy -> I sz, sy ** sx -> -I sz,
    sy ** sz -> I sx, sz ** sy -> -I sx,
    sz ** sx -> I sy, sx ** sz -> -I sy
  }
);

QACommutingRules[groups : {__List}] := (
  iSetNC[Flatten[groups]];
  Flatten @ Table[
    x ** y -> y ** x,
    {i, Length[groups]}, {j, i + 1, Length[groups]}, {x, groups[[j]]}, {y, groups[[i]]}
  ]
);


(* ::Section:: *)
(*Adjoints*)


Options[QADagger] = {"ComplexParameters" -> {}};

QADagger[expr_, modes : Except[_Rule | _RuleDelayed] : {}, OptionsPattern[]] := Module[
  {m, swap, cx, conj, wc},
  m = If[modes === {}, {}, iModes[modes]];
  cx = Alternatives @@ Flatten[{OptionValue["ComplexParameters"]}];
  conj = Function[k, (k /. Complex[re_, im_] :> Complex[re, -im]) /. s : cx :> Conjugate[s]];
  wc = If[m === $Failed, $Failed, iWordCoefficients[expr]];
  Which[
    m === $Failed, Message[QADagger::qamodes]; $Failed,
    wc === $Failed, Message[QADagger::qanonpoly]; $Failed,
    True,
      swap = Join[Rule @@@ m, Rule @@@ (Reverse /@ m)];
      NCExpand @ Total[KeyValueMap[conj[#2] iWord[Reverse[#1] /. swap] &, wc]]
  ]
];

Options[QAHermitianQ] = Join[Options[QADagger], Options[QANormalOrder]];
QAHermitianQ[expr_, modes_, rules_, opts : OptionsPattern[]] := Module[{d},
  d = QADagger[expr, modes, Sequence @@ FilterRules[{opts}, Options[QADagger]]];
  If[d === $Failed, $Failed,
    iZeroQ[QANormalOrder[expr - d, rules, Sequence @@ FilterRules[{opts}, Options[QANormalOrder]]]]
  ]
];


(* ::Section:: *)
(*Similarity transforms: Exp[lambda a] ** b ** Exp[-lambda a]*)


Options[QAExpandDualExponential] = {
  "Expansion" -> Automatic,
  "MaxOrder" -> 20,
  "ClosedForm" -> True,
  "MaxPeriod" -> 2,
  "MaxDimension" -> 20,
  "Assumptions" :> $Assumptions,
  "MaxIterations" -> 1000,
  "Method" -> Automatic
};


(* ::Subsection:: *)
(*Series method*)


(*
   Sums  b + lambda [a,b] + lambda^2/2! [a,[a,b]] + ...  with C_k = ad_a^k(b).
   Stops when C_n = 0, when C_n = mu C_(n-p) for some p <= "MaxPeriod" (then the tail
   is summed exactly), or after "MaxOrder" terms. If C_(k+p) = mu C_k for all k >= n-p, the tail is
     Sum_{r=n-p}^{n-1} C_r Sum_{j>=0} mu^j lambda^(r+p j)/(r+p j)!
   Returns {result, truncatedQ}.
*)
iAdjointSeries[a_, b_, rules_, lambda_, opts_List] := Module[
  {maxOrder, closedQ, maxPeriod, assum, normOpts, cs, c, n = 0, rec = None, trunc = False,
   mu, p, r, j, direct, tail},

  {maxOrder, closedQ, maxPeriod, assum} = OptionValue[
    QAExpandDualExponential, opts, {"MaxOrder", "ClosedForm", "MaxPeriod", "Assumptions"}];
  normOpts = FilterRules[opts, Options[QANormalOrder]];
  cs = {QANormalOrder[b, rules, Sequence @@ normOpts]};

  While[True,
    If[n >= maxOrder, trunc = True; Break[]];
    c = QACommutator[a, Last[cs], rules, Sequence @@ normOpts];
    n++;
    If[iZeroQ[c], Break[]];
    If[TrueQ[closedQ],
      Do[
        mu = iProportionality[c, cs[[-p]]];
        If[mu =!= None, rec = {p, mu}; Break[]],
        {p, Min[maxPeriod, n]}
      ];
      If[rec =!= None, Break[]]
    ];
    AppendTo[cs, c]
  ];

  If[rec === None,
    {NCExpand @ Total @ Table[lambda^k / k! cs[[k + 1]], {k, 0, Length[cs] - 1}], trunc},

    {p, mu} = rec;
    direct = Total @ Table[lambda^k / k! cs[[k + 1]], {k, 0, n - p - 1}];
    tail = Total @ Table[
      iSimplifyCoef[Sum[mu^j lambda^(r + p j) / (r + p j)!, {j, 0, Infinity}], assum] cs[[r + 1]],
      {r, n - p, n - 1}
    ];
    {NCExpand[direct + tail], False}
  ]
];


(* ::Subsection:: *)
(*Linear-closure method*)


(*
   Collect every word reachable from the words of bs by repeated [a, .]. If they span a space of
   at most "MaxDimension" words, d/d lambda W = M W, so W(lambda) = MatrixExp[lambda M] W(0).
   Returns the list of transformed operators, or $Failed.
*)
iLinearClosure[a_, bs_List, rules_, lambda_, opts_List] := Catch[Module[
  {maxDim, assum, normOpts, bcs, queue, rows = <||>, w, dc, basis, n, idx, mat, ex},

  {maxDim, assum} = OptionValue[QAExpandDualExponential, opts, {"MaxDimension", "Assumptions"}];
  normOpts = FilterRules[opts, Options[QANormalOrder]];

  bcs = iWordCoefficients[QANormalOrder[#, rules, Sequence @@ normOpts]] & /@ bs;
  If[MemberQ[bcs, $Failed], Throw[$Failed, iTag]];

  queue = DeleteDuplicates[Join @@ (Keys /@ bcs)];
  While[queue =!= {},
    w = First[queue]; queue = Rest[queue];
    If[!KeyExistsQ[rows, w],
      If[Length[rows] >= maxDim, Throw[$Failed, iTag]];
      dc = If[w === {}, <||>, iWordCoefficients[QACommutator[a, iWord[w], rules, Sequence @@ normOpts]]];
      If[dc === $Failed, Throw[$Failed, iTag]];
      rows[w] = dc;
      queue = Join[queue, Select[Keys[dc], !KeyExistsQ[rows, #] &]]
    ]
  ];

  basis = Keys[rows];
  n = Length[basis];
  idx = AssociationThread[basis -> Range[n]];
  mat = ConstantArray[0, {n, n}];
  KeyValueMap[
    Function[{w2, row}, KeyValueMap[(mat[[idx[w2], idx[#1]]] = #2) &, row]],
    rows
  ];
  ex = MatrixExp[lambda mat];

  Map[
    Function[bc,
      If[bc === <||>,
        0,
        NCExpand @ Total @ MapThread[
          iSimplifyCoef[#1, assum] iWord[#2] &,
          {Total[KeyValueMap[#2 ex[[idx[#1]]] &, bc]], basis}
        ]
      ]
    ],
    bcs
  ]
], iTag];


(* ::Subsection:: *)
(*Dispatch*)


iTransform[a_, b_, rules_, lambda_, caller_Symbol, opts_List] := Module[
  {method, maxOrder, maxDim, ser, lin},
  {method, maxOrder, maxDim} =
    OptionValue[QAExpandDualExponential, opts, {"Expansion", "MaxOrder", "MaxDimension"}];
  Switch[method,
    "LinearClosure",
      lin = iLinearClosure[a, {b}, rules, lambda, opts];
      If[lin === $Failed, Message[caller::closure, maxDim]; $Failed, First[lin]],
    "Series",
      ser = iAdjointSeries[a, b, rules, lambda, opts];
      If[ser[[2]], Message[caller::trunc, maxOrder]];
      ser[[1]],
    _,
      ser = iAdjointSeries[a, b, rules, lambda, opts];
      If[!ser[[2]],
        ser[[1]],
        lin = iLinearClosure[a, {b}, rules, lambda, opts];
        If[lin =!= $Failed, First[lin], Message[caller::trunc, maxOrder]; ser[[1]]]
      ]
  ]
];

QAExpandDualExponential[op_, a_, rules_, opts : OptionsPattern[]] :=
  iTransform[a, op, rules, 1, QAExpandDualExponential,
    FilterRules[Join[{opts}, Options[QAExpandDualExponential]], Options[QAExpandDualExponential]]];


(* ::Section:: *)
(*Dynamics*)


Options[QAPropagateOperator] = Append[Options[QAExpandDualExponential], "HBar" -> Global`\[HBar]];

(* X_H(t) = U^\[Dagger] X U with U = Exp[-I H t/\[HBar]], i.e. lambda = I t/\[HBar] and a = H. *)
QAPropagateOperator[op_, ham_, rules_, t_, opts : OptionsPattern[]] := (
  iCheckScalars[QAPropagateOperator, {t}];
  iTransform[ham, op, rules, I t / OptionValue["HBar"], QAPropagateOperator,
    FilterRules[Join[{opts}, Options[QAPropagateOperator]], Options[QAExpandDualExponential]]]
);

QAPropagateOperator[op_, ham_, rules_, t_, t0 : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
  QAPropagateOperator[op, ham, rules, t - t0, opts];

QAPropogateOperator[args___] := (Message[QAPropogateOperator::deprec]; QAPropagateOperator[args]);


Options[QAHeisenbergEOM] = Join[{"HBar" -> Global`\[HBar]}, Options[QANormalOrder]];
QAHeisenbergEOM[op_, ham_, rules_, opts : OptionsPattern[]] :=
  NCExpand[I / OptionValue["HBar"] *
    QACommutator[ham, op, rules, Sequence @@ FilterRules[{opts}, Options[QANormalOrder]]]];


Options[QAHeisenbergSolve] = Join[
  FilterRules[Options[QAExpandDualExponential], {"MaxDimension", "Assumptions", "MaxIterations", "Method"}],
  {"HBar" -> Global`\[HBar]}
];

QAHeisenbergSolve[ops_List, ham_, rules_, t_, opts : OptionsPattern[]] := Module[{res},
  iCheckScalars[QAHeisenbergSolve, {t}];
  res = iLinearClosure[ham, ops, rules, I t / OptionValue["HBar"],
    FilterRules[Join[{opts}, Options[QAHeisenbergSolve]], Options[QAExpandDualExponential]]];
  If[res === $Failed, Message[QAHeisenbergSolve::closure, OptionValue["MaxDimension"]]];
  res
];

QAHeisenbergSolve[op : Except[_List], ham_, rules_, t_, opts : OptionsPattern[]] :=
  With[{res = QAHeisenbergSolve[{op}, ham, rules, t, opts]}, If[res === $Failed, res, First[res]]];


(* ::Section:: *)
(*Baker-Campbell-Hausdorff*)


iNCPower[x_, 0] := 1;
iNCPower[x_, 1] := x;
iNCPower[x_, k_Integer] := NonCommutativeMultiply @@ ConstantArray[x, k];

(*
   With s a grading parameter, X(s) = Exp[s x] Exp[s y] - 1 = Sum_d X_d s^d and
   Z = log(1 + X) = Sum_k (-1)^(k+1) X^k / k, all truncated at s^n. Products are
   normal ordered as they are formed, which keeps the expressions small.
*)
QABCH[x_, y_, rules_, n : (_Integer?Positive) : 4, opts : OptionsPattern[QANormalOrder]] := Module[
  {norm, gmul, xs, pw, z},
  norm = QANormalOrder[#, rules, opts] &;
  gmul = Function[{u, v},
    Table[norm[Sum[u[[i + 1]] ** v[[d - i + 1]], {i, 0, d}]], {d, 0, n}]];
  xs = Table[
    If[d == 0, 0, norm[Sum[iNCPower[x, i] ** iNCPower[y, d - i] / (i! (d - i)!), {i, 0, d}]]],
    {d, 0, n}];
  pw = xs;
  z = xs;
  Do[
    pw = gmul[pw, xs];
    z = z + (-1)^(k + 1) / k pw,
    {k, 2, n}
  ];
  norm[Total[z]]
];


(* ::Section:: *)
(*Substitutions: displacement, squeezing, Bogoliubov*)


QAOperatorSubstitute[expr_, map_List, rules_, opts : OptionsPattern[QANormalOrder]] :=
  With[{wc = iWordCoefficients[expr]},
    If[wc === $Failed,
      Message[QAOperatorSubstitute::qanonpoly]; $Failed,
      QANormalOrder[Total[KeyValueMap[#2 iWord[#1 /. map] &, wc]], rules, opts]
    ]
  ];

QADisplace[expr_, modes_, alphas_, opts : OptionsPattern[QANormalOrder]] := Module[
  {m = iModes[modes], as = iList[alphas]},
  Which[
    m === $Failed, Message[QADisplace::qamodes]; $Failed,
    Length[as] =!= Length[m], Message[QADisplace::qalength]; $Failed,
    True,
      iCheckScalars[QADisplace, as];
      QAOperatorSubstitute[expr,
        Join[Thread[m[[All, 1]] -> m[[All, 1]] + as], Thread[m[[All, 2]] -> m[[All, 2]] + Conjugate[as]]],
        QABosonRules[m], opts]
  ]
];

QASqueeze[expr_, {a_, c_}, r_, theta : Except[_Rule | _RuleDelayed] : 0, opts : OptionsPattern[QANormalOrder]] := (
  iCheckScalars[QASqueeze, {r, theta}];
  QAOperatorSubstitute[expr,
    {a -> a Cosh[r] - c Exp[I theta] Sinh[r], c -> c Cosh[r] - a Exp[-I theta] Sinh[r]},
    QABosonRules[{a, c}], opts]
);


(* ::Section:: *)
(*Expectation values in bosonic states*)


(* {creator counts, annihilator counts} per mode for a normal-ordered word, or $Failed. *)
iBosonCounts[letters_, m_] := If[!SubsetQ[Flatten[m], letters],
  $Failed,
  {Count[letters, #] & /@ m[[All, 2]], Count[letters, #] & /@ m[[All, 1]]}
];

(* valueFn[ks, ls] gives the expectation value of Product_i c_i^k_i a_i^l_i. *)
iBosonExpectation[caller_Symbol, expr_, modes_, params_List, valueFn_] := Catch[Module[
  {m = iModes[modes], wc, parts},
  If[m === $Failed, Message[caller::qamodes]; Throw[$Failed, iTag]];
  If[!AllTrue[params, Length[#] == Length[m] &], Message[caller::qalength]; Throw[$Failed, iTag]];
  wc = iWordCoefficients[QANormalOrder[expr, QABosonRules[m]]];
  If[wc === $Failed, Message[caller::qanonpoly]; Throw[$Failed, iTag]];
  parts = KeyValueMap[
    With[{cnt = iBosonCounts[#1, m]}, If[cnt === $Failed, $Failed, #2 valueFn @@ cnt]] &,
    wc];
  If[MemberQ[parts, $Failed], Message[caller::ops]; Throw[$Failed, iTag]];
  Simplify[FunctionExpand[Total[parts]]]
], iTag];

QAFockExpectation[expr_, modes_, ns_] := With[{nl = iList[ns]},
  iBosonExpectation[QAFockExpectation, expr, modes, {nl},
    Function[{ks, ls}, Times @@ MapThread[If[#1 == #2, FactorialPower[#3, #1], 0] &, {ks, ls, nl}]]]
];

QAFockMatrixElement[expr_, modes_, ms_, ns_] := With[{ml = iList[ms], nl = iList[ns]},
  iBosonExpectation[QAFockMatrixElement, expr, modes, {ml, nl},
    Function[{ks, ls}, Times @@ MapThread[
      KroneckerDelta[#4, #3 - #2 + #1] Sqrt[FactorialPower[#3, #2] FactorialPower[#4, #1]] &,
      {ks, ls, nl, ml}]]]
];

QACoherentExpectation[expr_, modes_, alphas_] := With[{al = iList[alphas]},
  iBosonExpectation[QACoherentExpectation, expr, modes, {al},
    Function[{ks, ls}, Times @@ MapThread[Conjugate[#3]^#1 #3^#2 &, {ks, ls, al}]]]
];

QAThermalExpectation[expr_, modes_, nbars_] := With[{nb = iList[nbars]},
  iBosonExpectation[QAThermalExpectation, expr, modes, {nb},
    Function[{ks, ls}, Times @@ MapThread[If[#1 == #2, #1! #3^#1, 0] &, {ks, ls, nb}]]]
];


(* ::Section:: *)
(*Matrix representations*)


iIdentity[d_] := SparseArray[{i_, i_} -> 1, {d, d}];

iEmbed[mat_, i_, dims_] := If[Length[dims] == 1,
  mat,
  KroneckerProduct @@ Table[If[k == i, mat, iIdentity[dims[[k]]]], {k, Length[dims]}]
];

QAToMatrix[expr_, reps_List] := Catch[Module[{subs, dims, full, wc, used},
  subs = If[MatchQ[reps, {(_Rule | _RuleDelayed) ..}], {reps}, reps];
  If[!MatchQ[subs, {{(_Rule | _RuleDelayed) ..} ..}], Message[QAToMatrix::reps]; Throw[$Failed, iTag]];
  dims = Length[#[[1, 2]]] & /@ subs;
  full = Association @ Flatten @ MapIndexed[
    Function[{sub, pos}, (#[[1]] -> iEmbed[#[[2]], First[pos], dims]) & /@ sub],
    subs];
  wc = iWordCoefficients[expr];
  If[wc === $Failed, Message[QAToMatrix::qanonpoly]; Throw[$Failed, iTag]];
  used = Union @@ Keys[wc];
  If[!SubsetQ[Keys[full], used], Message[QAToMatrix::ops, Complement[used, Keys[full]]]; Throw[$Failed, iTag]];
  Total[KeyValueMap[#2 If[#1 === {}, iIdentity[Times @@ dims], Dot @@ Lookup[full, #1]] &, wc]] +
    0 iIdentity[Times @@ dims]
], iTag];

(* a|n> = Sqrt[n] |n-1>, basis |0>, |1>, ..., |d-1> *)
iBosonA[d_] := SparseArray[Table[{k, k + 1} -> Sqrt[k], {k, 1, d - 1}], {d, d}];

QABosonMatrices[modes_, nmax_] := Module[{m = iModes[modes], dims},
  If[m === $Failed,
    Message[QABosonMatrices::qamodes]; $Failed,
    iSetNC[Flatten[m]];
    dims = If[ListQ[nmax], nmax, ConstantArray[nmax, Length[m]]] + 1;
    Flatten @ Table[
      With[{am = iEmbed[iBosonA[dims[[i]]], i, dims]}, {m[[i, 1]] -> am, m[[i, 2]] -> Transpose[am]}],
      {i, Length[m]}]
  ]
];

(* Basis m = j, j-1, ..., -j.  J+|j m> = hbar Sqrt[j(j+1) - m(m+1)] |j m+1>. *)
iSpinLadder[j_,h_]:=Module[{d=2 j+1,ms=Range[j,-j,-1],jp,jm},
jp=SparseArray[Table[{i-1,i}->h Sqrt[j (j+1)-ms[[i]] (ms[[i]]+1)],{i,2,d}],{d,d}];
jm=SparseArray[Table[{i+1,i}->h Sqrt[j (j+1)-ms[[i]] (ms[[i]]-1)],{i,1,d-1}],{d,d}];
{jp,jm,SparseArray[h DiagonalMatrix[ms]]}];

iSpinQ[j_] := NumericQ[j] && TrueQ[j >= 0] && IntegerQ[2 j];

Options[QASpinMatrices] = {"HBar" -> Global`\[HBar]};
QASpinMatrices[{jx_, jy_, jz_}, j_?iSpinQ, OptionsPattern[]] :=
  Module[{jp, jm, jzm},
    iSetNC[{jx, jy, jz}];
    {jp, jm, jzm} = iSpinLadder[j, OptionValue["HBar"]];
    {jx -> (jp + jm) / 2, jy -> (jp - jm) / (2 I), jz -> jzm}
  ];

Options[QASpinLadderMatrices] = {"HBar" -> Global`\[HBar]};
QASpinLadderMatrices[{jp_, jm_, jz_}, j_?iSpinQ, OptionsPattern[]] :=
  (iSetNC[{jp, jm, jz}]; Thread[{jp, jm, jz} -> iSpinLadder[j, OptionValue["HBar"]]]);

QAPauliMatrices[{sx_, sy_, sz_}] :=
  (iSetNC[{sx, sy, sz}]; Thread[{sx, sy, sz} -> (SparseArray[PauliMatrix[#]] & /@ {1, 2, 3})]);


(* ::Section:: *)
(*Ordering prescriptions: symmetrization and Weyl ordering*)


(* Apply runFn to every maximal run of letters drawn from ops. runFn: run -> {{weight, newRun}, ...} *)
iReorder[expr_, ops_List, runFn_] :=
  NCExpand @ Total[iReorderTerm[#, ops, runFn] & /@ iTerms[NCExpand[expr]]];

iReorderTerm[t_, ops_, runFn_] := Module[{split = iSplitTerm[t], runs, choices},
  If[split === $Failed || split[[2]] === {},
    t,
    runs = SplitBy[split[[2]], MemberQ[ops, #] &];
    choices = If[MemberQ[ops, First[#]], runFn[#], {{1, #}}] & /@ runs;
    split[[1]] Total[(Times @@ #[[All, 1]]) iWord[Join @@ #[[All, 2]]] & /@ Tuples[choices]]
  ]
];

iSymmetrizeRun[a_, b_][run_] := {
  {1/2, Join[Select[run, # === a &], Select[run, # =!= a &]]},
  {1/2, Join[Select[run, # === b &], Select[run, # =!= b &]]}
};

iWeylRun[run_] := With[{perms = Permutations[run]}, {1 / Length[perms], #} & /@ perms];

QASymmetrize[expr_, a_, b_] := iReorder[expr, {a, b}, iSymmetrizeRun[a, b]];

QAWeylOrder[expr_, ops_List] := iReorder[expr, ops, iWeylRun];


(* ::Section:: *)
(*Angular momentum coupling*)


QACGTable[j1_?iSpinQ, j2_?iSpinQ, M_] /; NumericQ[M] && IntegerQ[j1 + j2 - M] && Abs[M] <= j1 + j2 :=
  Module[{m1s, js, data},
    m1s = Select[Range[j1, -j1, -1], Abs[M - #] <= j2 &];
    js = Range[j1 + j2, j1 + j2 - Length[m1s] + 1, -1];
    data = Table[ClebschGordan[{j1, m1}, {j2, M - m1}, {jj, M}], {m1, m1s}, {jj, js}];
    TableForm[data,
      TableHeadings -> {
        Row[{"(", #, ", ", M - #, ")"}] & /@ m1s,
        Row[{"J = ", #}] & /@ js
      }
    ]
  ];
QACGTable[___] := (Message[QACGTable::args]; $Failed);

QACGTables[j1_?iSpinQ, j2_?iSpinQ] :=
  Column[
    Labeled[QACGTable[j1, j2, #], Row[{"M = ", #}], Top] & /@ Range[j1 + j2, -(j1 + j2), -1],
    Spacings -> 2
  ];
QACGTables[___] := (Message[QACGTables::args]; $Failed);

QACoupledState[j1_?iSpinQ, j2_?iSpinQ, jj_, mm_] /;
    NumericQ[jj] && NumericQ[mm] && IntegerQ[j1 + j2 - jj] && Abs[j1 - j2] <= jj <= j1 + j2 &&
    IntegerQ[jj - mm] && Abs[mm] <= jj :=
  Association @ Select[
    Table[
      {m1, mm - m1} -> ClebschGordan[{j1, m1}, {j2, mm - m1}, {jj, mm}],
      {m1, Select[Range[j1, -j1, -1], Abs[mm - #] <= j2 &]}
    ],
    Last[#] =!= 0 &
  ];
QACoupledState[___] := (Message[QACoupledState::args]; $Failed);


(* ::Section:: *)
(*Display*)


iSetDaggerFormat[{x_, y_}] := (Format[y] := SuperDagger[x]);
iClearDaggerFormat[{_, y_}] := (FormatValues[y] = {});

QASetDaggerFormat[modes_] := With[{m = iModes[modes]},
  If[m === $Failed, Message[QASetDaggerFormat::qamodes]; $Failed, Scan[iSetDaggerFormat, m];]
];

QAClearDaggerFormat[modes_] := With[{m = iModes[modes]},
  If[m === $Failed, Message[QAClearDaggerFormat::qamodes]; $Failed, Scan[iClearDaggerFormat, m];]
];


End[];

SetAttributes[
  {QANormalOrder, QACommutator, QAAntiCommutator, QANestedCommutator, QACheckRules,
   QABosonRules, QAFermionRules, QACanonicalRules, QASpinRules, QASpinLadderRules, QAPauliRules,
   QACommutingRules, QADagger, QAHermitianQ,
   QAExpandDualExponential, QAPropagateOperator, QAPropogateOperator, QAHeisenbergEOM, QAHeisenbergSolve,
   QABCH, QAOperatorSubstitute, QADisplace, QASqueeze,
   QAFockExpectation, QAFockMatrixElement, QACoherentExpectation, QAThermalExpectation,
   QAToMatrix, QABosonMatrices, QASpinMatrices, QASpinLadderMatrices, QAPauliMatrices,
   QASymmetrize, QAWeylOrder, QACGTable, QACGTables, QACoupledState,
   QASetDaggerFormat, QAClearDaggerFormat},
  {Protected, ReadProtected}
];

EndPackage[];

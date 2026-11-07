(* ::Package:: *)

(* ::Title:: *)
(*QuantumAlgebra examples*)


(* ::Text:: *)
(*Open this file in Mathematica and evaluate it cell by cell, or run it all with Get.*)


(* ::Section:: *)
(*Setup*)


PacletDirectoryLoad[ParentDirectory[DirectoryName[$InputFileName /. "" :> NotebookFileName[]]]];
Needs["QuantumAlgebra`"]

(* Single lowercase letters are operators by default in NCAlgebra; these are scalars. *)
SetCommutative[t, m, n, r, g];


(* ::Section:: *)
(*Harmonic oscillator*)


bosons = QABosonRules[{a, c}];
QASetDaggerFormat[{a, c}];

QANormalOrder[a ** a ** c ** c, bosons]

(* Heisenberg picture: a(t) = E^(-I \[Omega] t) a *)
QAPropagateOperator[a, \[HBar] \[Omega] c ** a, bosons, t]

(* Position and momentum, solved as a linear system *)
canon = QACanonicalRules[{x, p}];
osc = p ** p/(2 m) + m \[Omega]^2 x ** x/2;
QAHeisenbergSolve[{x, p}, osc, canon, t, "Assumptions" -> \[HBar] > 0 && \[Omega] > 0 && m > 0]


(* ::Section:: *)
(*Spin precession*)


spin = QASpinRules[{jx, jy, jz}];
QAPropagateOperator[jx, \[Omega] jz, spin, t, "Assumptions" -> \[HBar] > 0 && \[Omega] > 0]

(* Check the algebra numerically for spin 3/2 *)
QAToMatrix[jx ** jx + jy ** jy + jz ** jz, QASpinMatrices[{jx, jy, jz}, 3/2]] // MatrixForm


(* ::Section:: *)
(*Beam splitter (two modes)*)


modes = {{a1, c1}, {a2, c2}};
bs = QABosonRules[modes];
QASetDaggerFormat[modes];

QAPropagateOperator[a1, \[HBar] g (c1 ** a2 + c2 ** a1), bs, t, "Assumptions" -> g > 0]

(* Hong-Ou-Mandel: at a 50/50 splitter, <1,1| n1 n2 |1,1> vanishes in the output *)
out = Thread[{a1, a2} -> QAHeisenbergSolve[{a1, a2}, \[HBar] g (c1 ** a2 + c2 ** a1), bs, \[Pi]/(4 g)]];
outDag = Thread[{c1, c2} -> (QADagger[#, modes] & /@ ({a1, a2} /. out))];
QAFockExpectation[
  QAOperatorSubstitute[c1 ** a1 ** c2 ** a2, Join[out, outDag], bs],
  modes, {1, 1}]


(* ::Section:: *)
(*Jaynes-Cummings model*)


pauli = QAPauliRules[{sx, sy, sz}];
jc = Join[bosons, pauli, QACommutingRules[{{a, c}, {sx, sy, sz}}]];
QACheckRules[jc]["Confluent"]

sp = (sx + I sy)/2; sm = (sx - I sy)/2;
hJC = \[HBar] \[Omega] c ** a + \[HBar] \[Omega] sz/2 + \[HBar] g (c ** sm + a ** sp);

(* The excitation number commutes with the Hamiltonian *)
QACommutator[hJC, c ** a + sz/2, jc]

QAHermitianQ[hJC, {a, c}, jc]


(* ::Section:: *)
(*Coherent and squeezed states*)


SetCommutative[\[Alpha]];

(* Displacement operator from BCH: e^(\[Alpha] c) e^(-\[Alpha]* a) = e^(\[Alpha] c - \[Alpha]* a + |\[Alpha]|^2/2) *)
QABCH[\[Alpha] c, -Conjugate[\[Alpha]] a, bosons]

(* Photon-number variance of a coherent state equals the mean *)
QACoherentExpectation[c ** a ** c ** a, {a, c}, \[Alpha]] - QACoherentExpectation[c ** a, {a, c}, \[Alpha]]^2 // Simplify

(* Quadrature variance of squeezed vacuum *)
xq = (a + c)/Sqrt[2];
QAFockExpectation[QASqueeze[xq ** xq, {a, c}, r], {a, c}, 0] // Simplify


(* ::Section:: *)
(*Ordering prescriptions*)


QAWeylOrder[x ** x ** p, {x, p}]
QANormalOrder[%, canon]


(* ::Section:: *)
(*Angular momentum coupling*)


QACGTables[1, 1/2]
QACoupledState[1, 1, 0, 0]

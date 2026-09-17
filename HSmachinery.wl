(* ::Package:: *)

(* ================================================================ *)
(*  Package: HSmachinery.m                                          *)
(*  Description: Tools for working with HS generating systems       *)
(*  Authors: Abrogolat, Qwen3.7, DeepSeek AI                        *)
(*  Date: 2026-09-16                                                *)
(*  Version: 1.1.0                                                  *)
(* ================================================================ *)

BeginPackage["HSmachinery`"]

(* --- Public API Usage Messages --- *)
spinors::usage = "spinors[vars] declares vars as spinor variables for use in convolutions and integrations.";

conv::usage = "conv[a, b] is the antisymmetric Sp(2) convolution of two spinors a and b.";

comp::usage = "comp[a, i] represents the i-th component of spinor a (i = 1, 2).";

integratePair::usage = 
"integratePair[{polynom, exponent}, {u, v}] integrates polynom * Exp[I * exponent] \
over the spinor pair u and v.";

reduceToSpinorInvariants::usage = 
"reduceToSpinorInvariants[poly] reduces a polynomial in spinor components to \
Sp(2)-invariant form using conv expressions.";

starProduct::usage = 
"starProduct[{poly1, exp1}, {poly2, exp2}, {z, y}, {zb, yb}] computes the \
Sp(2)-covariant star product of two phase-space symbols. Each symbol is given as \
{polynomial, exponent}, where both are expressions in spinor components. The \
spinors z, y, zb, yb must already be declared via spinors[]. The result is \
{integratedPolynomial, overallPhase} after performing Gaussian integration over \
auxiliary spinor pairs.";

starProductHol::usage = 
"starProduct[{poly1, exp1}, {poly2, exp2}, {z, y}] computes the holomorphic\
Sp(2)-covariant star product of two phase-space symbols. Each symbol is given as \
{polynomial, exponent}, where both are expressions in spinor components. The \
spinors z, y must already be declared via spinors[]. The result is \
{integratedPolynomial, overallPhase} after performing Gaussian integration over \
auxiliary spinor pairs.";

homotopy::usage =
  "Homotopy[{poly, exp}, {z, \[Theta]}, {\[Tau], q}, degree] computes the shifted \
homotopy of polynomial poly times exponential expression exp in the \
spinor variables z and \[Theta], with homotopy parameter \[Tau] and \
spinor shift q. degree must be 1 (linear in \[Theta]) or 2 (quadratic in \[Theta]). z, \[Theta] and q \
must belong to $Spinors.";

extractExponent::usage = "blah-blah";

Begin["`Private`"]

(* --- Initialization --- *)
$Spinors = {};

(* --- Public Function Implementations --- *)
spinors[vars__] := ($Spinors = Union[$Spinors, {vars}];)

(* Dummy symbol used to mask conv objects during scalar detection *)
$scalarDummy = Unique["scalarDummy"];

(* ScalarQ[expr] returns True if expr contains no spinor variables
   outside of conv[...] objects.  Since all spinors appear only through
   conv, we can safely replace every conv[__] by a dummy and then check
   for the spinor variables. *)
ScalarQ[expr_] := FreeQ[
    expr /. conv[__] -> $scalarDummy,
    Alternatives @@ $Spinors
];

(* General spinor test *)
SpinorQ[{xs__}] := VectorQ[{xs}, SpinorQ];   (* non-empty list: all elements *)
SpinorQ[expr_]  := !ScalarQ[expr];           (* everything else *)
SpinorQ[{}]     := False;                    (* {} is not spinorial *)

(* Atomic: x is literally one of the declared spinor symbols *)
AtomicSpinorQ[x_]      := MemberQ[$Spinors, x];
AtomicSpinorQ[{xs__}]  := VectorQ[{xs}, AtomicSpinorQ];
AtomicSpinorQ[_]       := False;   (* e.g. 2 z, z + y, conv[z,y], or {} *)

(* ---------- conv: all rules ----------- *)
(* antisymmetric, bilinear *)
conv[a_, a_] := 0
conv[a_ + b_, c_] := conv[a, c] + conv[b, c]
conv[a_, b_ + c_] := conv[a, b] + conv[a, c]

(* Coefficient extraction using ScalarQ *)
conv[c_ a_, b_] := c conv[a, b] /; ScalarQ[c]
conv[a_, c_ b_] := c conv[a, b] /; ScalarQ[c]

(* Antisymmetry (placed after coefficient extraction for efficiency) *)
conv[a_, b_] := -conv[b, a] /; !OrderedQ[{a, b}]
conv[a_, 0] := 0
conv[0, a_] := 0

(* ---------- comp: linearity ---------- *)
comp[a_ + b_, i_] := comp[a, i] + comp[b, i]
comp[c_ a_, i_] := c comp[a, i] /; ScalarQ[c]
comp[0, i_] := 0

(* ---------- Epsilon tensor ---------- *)
eps[1, 2] = 1;
eps[2, 1] = -1;
eps[1, 1] = 0;
eps[2, 2] = 0;

(* ---------- Helper: Extract T, A, B, C from exponent ---------- *)
extractExponent[exp_, u_, v_] := Module[
  {expanded, terms, T, A, B, C, a, b, c},

  expanded = Expand[exp];

  (* Extract all terms that are a scalar coefficient times a single conv *)
  terms = Cases[expanded,
    c_. * conv[a_, b_] /; FreeQ[c, conv] :> {c, a, b},
    {1}
  ];

  T = 0;
  A = 0;
  B = 0;

  Do[
    {c, a, b} = term;
    Which[
      a === u && b === v, T += c,
      a === v && b === u, T -= c,
      a === u && b =!= v, A += c * b,
      a =!= v && b === u, A -= c * a,
      a === v && b =!= u, B -= c * b,
      a =!= u && b === v, B += c * a,
      True, Null
    ],
    {term, terms}
  ];

  (* Whatever is left over goes into C *)
  C = Expand[expanded - T conv[u, v] - conv[u, A] - conv[B, v]];

  {T, A, B, C}
];


(* ---------- Helper: Convert conv to components ---------- *)
convToComponents[expr_] := expr /. conv[a_, b_] :> comp[a, 1] comp[b, 2] - comp[a, 2] comp[b, 1]


(* ---------- Helper: Generate all monomials of degree k from basis ---------- *)
generateMonomials[basis_List, k_Integer] := Module[{n, expLists, monomials},
  n = Length[basis];

  If[k == 0, Return[{1}]];
  If[n == 0, Return[{}]];

  (* Generate all exponent lists that sum to k *)
  expLists = Flatten[
    Table[
      Flatten[Permutations[PadRight[part, n]]],
      {part, Join[{{k}}, IntegerPartitions[k]]}
    ],
    1
  ];

  (* Remove duplicates and generate monomials *)
  expLists = DeleteDuplicates[expLists];
  monomials = Times @@ (basis^#) & /@ expLists;

  Union[monomials]
]



reduceToSpinorInvariants[poly_] := Module[
  {p, compVars, spinorsInPoly, terms, homogeneousParts, result = 0,
   d, k, pairs, basis, monomials, target, basisMonomialsComp,
   allMonomials, eqns, unknowns, mat, vec, solution, sol,
   monomialsOf, getCoeff},

  p = Expand[poly];
  If[p === 0, Return[0]];

  compVars = Union[Cases[p, comp[_, _], Infinity]];
  spinorsInPoly = Union[Cases[p, comp[X_, _] :> X, Infinity]];

  If[Length[spinorsInPoly] == 0, Return[p]];

  terms = If[Head[p] === Plus, List @@ p, {p}];
  homogeneousParts = GroupBy[terms, Total[Exponent[#, compVars]] &];

  (* Robust extraction of monomials in compVars *)
  monomialsOf[expr_] := Module[{rules},
    rules = CoefficientRules[Expand[expr], compVars];
    Union[Times @@@ (compVars^# & /@ rules[[All, 1]])]
  ];

  getCoeff[expr_, mon_] := If[mon === 1,
    Coefficient[expr, compVars, ConstantArray[0, Length[compVars]]],
    Coefficient[expr, mon]
  ];

  Do[
    d = deg;
    k = d/2;

    If[! IntegerQ[k] || k < 0,
      result = result + Total[homogeneousParts[deg]];
      Continue[];
    ];

    pairs = Select[Subsets[spinorsInPoly, {2}], OrderedQ];

    If[Length[pairs] == 0 && d > 0,
      result = result + Total[homogeneousParts[deg]];
      Continue[];
    ];

    basis = conv @@@ pairs;

    If[k == 0,
      basisMonomialsComp = {1};
      monomials = {1};
      ,
      monomials = Union[Times @@@ Tuples[basis, k]];
      basisMonomialsComp =
        Expand[# /. conv[a_, b_] :>
          comp[a, 1] comp[b, 2] - comp[a, 2] comp[b, 1]] & /@ monomials;
    ];

    target = Expand[Total[homogeneousParts[deg]]];

    allMonomials = Union[
      monomialsOf[target],
      Flatten[monomialsOf /@ basisMonomialsComp]
    ];

    vec = getCoeff[target, #] & /@ allMonomials;
    mat = Table[
      getCoeff[basisMonomialsComp[[j]], allMonomials[[i]]],
      {i, Length[allMonomials]}, {j, Length[monomials]}
    ];

    If[Length[monomials] == 0,
      solution = {};
      ,
      unknowns = Array[c, Length[monomials]];
      eqns = Thread[mat . unknowns == vec];
      sol = Solve[eqns, unknowns];

      If[sol === {},
        result = result + target;
        Continue[];
      ];

      solution = unknowns /. First[sol];

      (* Pick one particular solution if the system is underdetermined *)
      solution = solution /. c[_] :> 0;
    ];

    result = result + (solution . monomials);

    , {deg, Keys[homogeneousParts]}];

  Expand[result]
];



(* ---------- Core: Gaussian integration over u, v ---------- *)
integrateGaussian[polynom_, u_, v_, T_] := Module[
  {expanded, terms, result = 0, term, monomials, mon, U, V, external},

  expanded = Expand[polynom];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];

  Do[
    term = terms[[k]];

    (* Replace conv with explicit component form for integration *)
    term = term /. {
      conv[u, v] :> u[1] v[2] - u[2] v[1],
      conv[u, X_] :> u[1] comp[X, 2] - u[2] comp[X, 1],
      conv[v, Y_] :> v[1] comp[Y, 2] - v[2] comp[Y, 1]
    };
    term = Expand[term];
    monomials = If[Head[term] === Plus, List @@ term, {term}];

    Do[
      mon = monomials[[m]];
      If[mon === 0, Continue[]];

      U = Table[Exponent[mon, u[i]], {i, 1, 2}];
      V = Table[Exponent[mon, v[i]], {i, 1, 2}];

      (* Non-zero only if U1 == V2 and U2 == V1 *)
      If[U[[1]] == V[[2]] && U[[2]] == V[[1]],
        external = mon /. {u[1] -> 1, u[2] -> 1, v[1] -> 1, v[2] -> 1};
        result = result + (1/T^2) * (I/T)^(U[[1]] + U[[2]]) * (-1)^U[[2]] * U[[1]]! * U[[2]]! * external;
      ];
    , {m, 1, Length[monomials]}];
  , {k, 1, Length[terms]}];

  (* Rigorously reduce the component polynomial back to Sp(2) invariants *)
  result = reduceToSpinorInvariants[result];

  result
];

(* ---------- Public API: integratePair ---------- *)
integratePair[{polynom_, exponent_}, {u_Symbol, v_Symbol}] := Module[
  {tVal, aVal, bVal, cVal, shiftedPoly, intResult, expFactor},

  {tVal, aVal, bVal, cVal} = extractExponent[exponent, u, v];

  If[tVal === 0,
    Message[integratePair::zeroT];
    Return[$Failed];
  ];

  (* Shift: u -> u - bVal/tVal, v -> v - aVal/tVal *)
  shiftedPoly = polynom /. {u -> u - bVal/tVal, v -> v - aVal/tVal};

  intResult = integrateGaussian[shiftedPoly, u, v, tVal];
  expFactor = cVal - conv[bVal, aVal]/tVal;

  (* Use Map to apply Expand to each element, avoiding direct assignment issues *)
  Map[Expand, {intResult, expFactor}]
];

integratePair::badvars = "Integration variables must be a list of two unevaluated symbols, e.g., {u, v}. Use Clear[u, v] if needed.";
integratePair::zeroT = "The coefficient T of conv[u, v] in the exponent is zero. Gaussian integration requires T != 0.";


(* ---------- Public API: starProduct ---------- *)
(* Declare temporary spinors locally without affecting global $Spinors *)
With[{tempSpinors = {s, t, sb, tb}},
  starProduct[
    {poly1_, exp1_},
    {poly2_, exp2_},
    {z_, y_},
    {zb_, yb_}
  ] := Module[{
      oldSpinors = $Spinors,
      shiftedPoly1, shiftedPoly2,
      shiftedExp1, shiftedExp2,
      combinedPoly, combinedExp,
      innerResult, outerResult
    },
    
    (* Temporarily extend spinor list *)
    $Spinors = Join[$Spinors, tempSpinors];
    
    (* Apply shifts to first factor: +s, +sb *)
    shiftedPoly1 = poly1 /. {
      z -> z + s,
      y -> y + s,
      zb -> zb + sb,
      yb -> yb + sb
    };
    shiftedExp1 = exp1 /. {
      z -> z + s,
      y -> y + s,
      zb -> zb + sb,
      yb -> yb + sb
    };
    
    (* Apply shifts to second factor: -t, +t, -tb, +tb *)
    shiftedPoly2 = poly2 /. {
      z -> z - t,
      y -> y + t,
      zb -> zb - tb,
      yb -> yb + tb
    };
    shiftedExp2 = exp2 /. {
      z -> z - t,
      y -> y + t,
      zb -> zb - tb,
      yb -> yb + tb
    };
    
    (* Combine polynomial and exponent *)
    combinedPoly = Expand[shiftedPoly1 * shiftedPoly2];
    combinedExp = Expand[
      shiftedExp1 + shiftedExp2 + conv[s, t] + conv[sb, tb]
    ];
    
    (* Perform inner integration over {sb, tb} *)
    innerResult = integratePair[{combinedPoly, combinedExp}, {sb, tb}];
    If[innerResult === $Failed,
      $Spinors = oldSpinors;
      Return[$Failed]
    ];
    
    (* Now integrate over {s, t} *)
    outerResult = integratePair[innerResult, {s, t}];
    
    (* Restore original spinor list *)
    $Spinors = oldSpinors;
    
    (* Return final result *)
    outerResult
  ] /; VectorQ[{z, y, zb, yb}, MemberQ[$Spinors, #] &]
];

(* Fallback message if undeclared spinors are passed *)
starProduct[___] := (
  Message[starProduct::invspinors, "All spinors in {{z,y},{zb,yb}} must be declared via spinors[]."];
  $Failed
);
starProduct::invspinors = "`1`";

(* ---------- Public API: starProductHol ---------- *)
(* Declare temporary spinors locally without affecting global $Spinors *)
With[{tempSpinors = {s, t}},
  starProductHol[
    {poly1_, exp1_},
    {poly2_, exp2_},
    {z_, y_}
  ] := Module[{
      oldSpinors = $Spinors,
      shiftedPoly1, shiftedPoly2,
      shiftedExp1, shiftedExp2,
      combinedPoly, combinedExp,
      innerResult, outerResult
    },
    
    (* Temporarily extend spinor list *)
    $Spinors = Join[$Spinors, tempSpinors];
    
    (* Apply shifts to first factor: +s, +sb *)
    shiftedPoly1 = poly1 /. {
      z -> z + s,
      y -> y + s
    };
    shiftedExp1 = exp1 /. {
      z -> z + s,
      y -> y + s
    };
    
    (* Apply shifts to second factor: -t, +t, -tb, +tb *)
    shiftedPoly2 = poly2 /. {
      z -> z - t,
      y -> y + t
    };
    shiftedExp2 = exp2 /. {
      z -> z - t,
      y -> y + t
    };
    
    (* Combine polynomial and exponent *)
    combinedPoly = Expand[shiftedPoly1 * shiftedPoly2];
    combinedExp = Expand[
      shiftedExp1 + shiftedExp2 + conv[s, t]
    ];
    
    (* Now integrate over {s, t} *)
    outerResult = integratePair[{combinedPoly, combinedExp}, {s, t}];
    
    (* Restore original spinor list *)
    $Spinors = oldSpinors;
    
    (* Return final result *)
    outerResult
  ] /; VectorQ[{z, y}, MemberQ[$Spinors, #] &]
];

(* Fallback message if undeclared spinors are passed *)
starProductHol[___] := (
  Message[starProduct::invspinors, "All spinors in {{z,y}} must be declared via spinors[]."];
  $Failed
);
starProduct::invspinors = "`1`";



(* ---------- Public API: Homotopy ---------- *)
homotopy[{poly_, exp_}, {z_, \[Theta]_}, {\[Tau]_, q_}, degree_Integer] /;
  SpinorQ[q]&&AtomicSpinorQ[{z, \[Theta]}] :=
  Module[{result},
    Which[
      degree === 1,
        {(poly/.{z-> \[Tau] z - (1-\[Tau])q, \[Theta] -> z + q}), (exp/.{z-> \[Tau] z - (1-\[Tau])q})},

      degree === 2,
        {2 \[Tau] comp[\[Theta], z + q] *(poly/.{z-> \[Tau] z - (1-\[Tau])q}), (exp/.{z-> \[Tau] z - (1-\[Tau])q})},

      True,
        Message[Homotopy::badDegree, degree];
        result = $Failed
    ];
    result
  ];

Homotopy[args___] := (
  Message[Homotopy::badArgs, HoldForm[Homotopy[args]]];
  $Failed
);
 
Homotopy::badDegree =
  "Homotopy degree `1` is not supported. Only degree = 1 or degree = 2 are allowed.";
  
Homotopy::badArgs =
  "Arguments `1` do not match the expected signature \
Homotopy[{poly, exp}, {z, \[Theta]}, {\[Tau], p}, degree].";


End[]
EndPackage[]


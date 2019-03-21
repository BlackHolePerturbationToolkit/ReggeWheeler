BeginPackage["ReggeWheeler`NumericalIntegration`"];


Begin["`Private`"];

SetAttributes[PsiIn, {NumericFunction}];

PsiIn[s_, l_, \[Omega]_, rmin_, rmax_] :=
 Module[{},
];

Derivative[1][PsiIn[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][r_?InexactNumberQ] :=
 Module[{},
];

Derivative[n_Integer?Positive][PsiIn[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][r0_?InexactNumberQ] :=
 Module[{d2R, Rderivs, R, r, i},
  d2R = Derivative[1][R][r] + R[r];

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> PsiIn[s, l, \[Omega]]'[r0],
    R[r] -> PsiIn[s, l, \[Omega]][r0], r -> r0}
];


SetAttributes[PsiUp, {NumericFunction}];

PsiUp[s_, l_, \[Omega]_, rmin_, rmax_] :=
 Module[{},
];

Derivative[1][PsiUp[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][r_?InexactNumberQ] :=
 Module[{},
];

Derivative[n_Integer?Positive][PsiUp[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][r0_?InexactNumberQ] :=
 Module[{d2R, Rderivs, R, r, i},
  d2R = Derivative[1][R][r] + R[r];

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> PsiUp[s, l, \[Omega]]'[r0],
    R[r] -> PsiUp[s, l, \[Omega]][r0], r -> r0}
];

End[]
EndPackage[];

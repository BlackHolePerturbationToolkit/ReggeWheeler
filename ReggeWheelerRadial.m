(* ::Package:: *)

BeginPackage["ReggeWheeler`ReggeWheelerRadial`",
  {
    "ReggeWheeler`NumericalIntegration`",
    "ReggeWheeler`MST`RenormalizedAngularMomentum`",
    "ReggeWheeler`MST`MST`",
    "SpinWeightedSpheroidalHarmonics`"
  }];

ReggeWheelerRadial::usage = "ReggeWheelerRadial[s, l, \[Omega]] computes solutions to the Regge Wheeler equation."
ReggeWheelerRadialFunction::usage = "ReggeWheelerRadialFunction[s, l, \[Omega], assoc] an object representing solutions to the Regge Wheeler equation."

Begin["`Private`"];

Options[ReggeWheelerRadial] = {Method -> {"NumericalIntegration", "rmin" -> 4, "rmax" -> 20}, "BoundaryConditions" -> {"In","Up"}};

ReggeWheelerRadial[s_Integer, l_Integer, \[Omega]_, OptionsPattern[]] :=
 Module[{assocIn, assocUp, \[Lambda], \[Nu], solFuncs, method, norms, m = 0, a=0},
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  Switch[OptionValue[Method],
    "MST"|{"MST", ___},
      \[Nu] = RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method->"RenormalizedAngularMomentum"/.OptionValue[Method][[2;;]] ];
      method = {"MST", "RenormalizedAngularMomentum" -> \[Nu]};
      norms = ReggeWheeler`MST`MST`Private`Amplitudes[s,l,m,a,2\[Omega],\[Nu],\[Lambda]];
      solFuncs =
       <|"In" -> ReggeWheeler`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In"]["Transmission"]],
         "Up" -> ReggeWheeler`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up"]["Transmission"]]|>;
      norms = norms/norms[[All, "Transmission"]];,
    {"NumericalIntegration", "rmin" -> _, "rmax" -> _},
      method = Association[OptionValue[Method][[2;;]]];
      solFuncs = 
       <|"In" -> ReggeWheeler`NumericalIntegration`Private`PsiIn[s, l, \[Omega], method["rmin"], method["rmax"]]["Psi"],
         "Up" -> ReggeWheeler`NumericalIntegration`Private`PsiUp[s, l, \[Omega], method["rmin"], method["rmax"]]["Psi"]|>;
  ];

  assocIn = Association[
    "s" -> s,
    "l" -> l,
    "m" -> m,
    "\[Omega]" -> \[Omega],
    "Method" -> method,
    "BoundaryConditions" -> OptionValue["BoundaryConditions"],
    "Eigenvalue" -> \[Lambda],
    "Amplitudes" -> norms["In"],
    "RadialFunction" -> solFuncs["In"]
    ];

  assocUp = Association[
    "s" -> s,
    "l" -> l,
    "m" -> m,
    "\[Omega]" -> \[Omega],
    "Method" -> method,
    "BoundaryConditions" -> "Up",
    "Eigenvalue" -> \[Lambda],
    "Amplitudes" -> norms["Up"],
    "RadialFunction" -> solFuncs["Up"]
    ];

  <|"In" -> ReggeWheelerRadialFunction[s, l, \[Omega], assocIn],
    "Up" -> ReggeWheelerRadialFunction[s, l, \[Omega], assocUp]|>
];

Format[ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]] := 
  "ReggeWheelerRadialFunction["<>ToString[s]<>","<>ToString[l]<>","<>ToString[\[Omega]]<>",<<>>]";

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][y_String] /; !MemberQ[{"RadialFunction"}, y] :=
  assoc[y];

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][r_?NumericQ] :=
  assoc["RadialFunction"][r];

Derivative[n_][ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]][r_?NumericQ] :=
  assoc["RadialFunction"]'[r];


End[]
EndPackage[];

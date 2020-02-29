(* ::Package:: *)

BeginPackage["ReggeWheeler`ReggeWheelerRadial`",
  {
    "ReggeWheeler`NumericalIntegration`",
    "ReggeWheeler`MST`RenormalizedAngularMomentum`",
    "ReggeWheeler`MST`MST`"
  }];

ReggeWheelerRadial::usage = "ReggeWheelerRadial[s, l, \[Omega]] computes solutions to the Regge Wheeler equation."
ReggeWheelerRadialFunction::usage = "ReggeWheelerRadialFunction[s, l, \[Omega], assoc] an object representing solutions to the Regge Wheeler equation."

Begin["`Private`"];

Options[ReggeWheelerRadial] = {Method -> {"NumericalIntegration", "rmin" -> 4, "rmax" -> 20}, "BoundaryConditions" -> {"In","Up"}};

ReggeWheelerRadial[s_Integer, l_Integer, \[Omega]_, OptionsPattern[]] :=
 Module[{assoc, \[Lambda], \[Nu], solFuncs, method, norms, m = 0, a=0},
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  Switch[OptionValue[Method],
    "MST"|{"MST", ___},
      \[Nu] = RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method->"RenormalizedAngularMomentum"/.OptionValue[Method][[2;;]] ];
      method = {"MST", "RenormalizedAngularMomentum" -> \[Nu]};
      norms = ReggeWheeler`MST`MST`Private`Amplitudes[s,l,m,a,2\[Omega],\[Nu],\[Lambda]];
      solFuncs = OptionValue["BoundaryConditions"] /.
        {"In" -> ReggeWheeler`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In"]["Transmission"]],
         "Up" -> ReggeWheeler`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up"]["Transmission"]]};
      norms = norms/norms[[All, "Transmission"]];,
    {"NumericalIntegration", "rmin" -> _, "rmax" -> _},
      method = Association[OptionValue[Method][[2;;]]];
      solFuncs = 
        {ReggeWheeler`NumericalIntegration`Private`PsiIn[s, l, \[Omega], method["rmin"], method["rmax"]],
         ReggeWheeler`NumericalIntegration`Private`PsiUp[s, l, \[Omega], method["rmin"], method["rmax"]]};
  ];

  assoc = Association[
    "s" -> s,
    "l" -> l,
    "m" -> m,
    "\[Omega]" -> \[Omega],
    "Method" -> method,
    "BoundaryConditions" -> OptionValue["BoundaryConditions"],
    "Eigenvalue" -> \[Lambda],
    "Amplitudes" -> norms,
    "SolutionFunctions" -> solFuncs
    ];

  ReggeWheelerRadialFunction[s, l, \[Omega], assoc]
];

Format[ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]] := 
  "ReggeWheelerRadialFunction["<>ToString[s]<>","<>ToString[l]<>","<>ToString[\[Omega]]<>",<<>>]";

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][y:("In"|"Up")] :=
 Module[{assocNew=assoc},
  assocNew["SolutionFunctions"] = First[Pick[assoc["SolutionFunctions"], assoc["BoundaryConditions"], y]];
  assocNew["BoundaryConditions"] = y;
  ReggeWheelerRadialFunction[s, l, \[Omega], assocNew]
];

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][y_String] /; !MemberQ[{"SolutionFunctions"},y]:= assoc[y];

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][r_?NumericQ] := Module[{},
  If[
    Head[assoc["BoundaryConditions"]] === List,
    Return[Association[MapThread[#1 -> #2["Psi"][r] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]], 
    Return[assoc["SolutionFunctions"]["Psi"][r]]
  ];  
];

Derivative[n_][ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]][r_?NumericQ] := Module[{},
  If[
    Head[assoc["BoundaryConditions"]] === List,
    Return[Association[MapThread[#1 -> #2["dPsidr"][r] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]], 
    Return[assoc["SolutionFunctions"]["dPsidr"][r]]
  ];  
];


End[]
EndPackage[];

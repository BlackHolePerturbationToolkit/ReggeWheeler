BeginPackage["ReggeWheeler`ReggeWheelerRadial`",
  {
    "ReggeWheeler`NumericalIntegration`"
  }];

ReggeWheelerRadial::usage = "ReggeWheelerRadial[s, l, \[Omega]] computes solutions to the Regge Wheeler equation."
ReggeWheelerRadialFunction::usage = "ReggeWheelerRadialFunction[s, l, \[Omega], assoc] an object representing solutions to the Regge Wheeler equation."

Begin["`Private`"];

Options[ReggeWheelerRadial] = {Method -> {"NumericalIntegration", "rmin" -> 4, "rmax" -> 20}, "BoundaryConditions" -> {"In", "Up"}};

ReggeWheelerRadial[s_Integer, l_Integer, \[Omega]_, OptionsPattern[]] :=
 Module[{assoc, rmin, rmax, solFuncs, method},
  Switch[OptionValue[Method],
    "MST"|{"MST", ___},
      method = {"MST"};
      solFuncs = $Failed,
    {"NumericalIntegration", "rmin" -> _, "rmax" -> _},
      method = OptionValue[Method];
      solFuncs = OptionValue["BoundaryConditions"] /.
        {"In" -> ReggeWheeler`NumericalIntegration`Private`PsiIn[s, l, \[Omega], rmin, rmax],
         "Up" -> ReggeWheeler`NumericalIntegration`Private`PsiUp[s, l, \[Omega], rmin, rmax]};
  ];

  assoc = Association[
    "Method" -> method,
    "BoundaryConditions" -> OptionValue["BoundaryConditions"],
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
    Return[Association[MapThread[#1 -> #2[r] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]], 
    Return[assoc["SolutionFunctions"][r]]
  ];  
];

Derivative[n_][ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]][r_?NumericQ] := Module[{},
  If[
    Head[assoc["BoundaryConditions"]] === List,
    Return[Association[MapThread[#1 -> Derivative[n][#2][r] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]], 
    Return[Derivative[n][assoc["SolutionFunctions"]][r]]
  ];  
];


End[]
EndPackage[];

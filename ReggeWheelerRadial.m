(* ::Package:: *)

(* ::Title:: *)
(*ReggeWheeler package*)


(* ::Section::Closed:: *)
(*Create Package*)


BeginPackage["ReggeWheeler`ReggeWheelerRadial`",
  {
    "ReggeWheeler`NumericalIntegration`",
    "ReggeWheeler`MST`RenormalizedAngularMomentum`",
    "ReggeWheeler`MST`MST`",
    "SpinWeightedSpheroidalHarmonics`"
  }];


(* Usage messages *)
ReggeWheelerRadial::usage = "ReggeWheelerRadial[s, l, \[Omega]] computes homogeneous solutions to the Regge Wheeler equation."
ReggeWheelerRadialFunction::usage = "ReggeWheelerRadialFunction[s, l, \[Omega], assoc] is an object representing a homogeneous solution to the Regge Wheeler equation."


(* Error messages *)
ReggeWheelerRadial::optx = "Unknown options in `1`";


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ReggeWheelerRadial*)


(* ::Subsection::Closed:: *)
(*Numerical Integration Method*)


Options[ReggeWheelerRadialNumericalIntegration] = {
  "rmin" -> None,
  "rmax" -> None,
  "BoundaryConditions" -> None};


ReggeWheelerRadialNumericalIntegration[s_Integer, l_Integer, \[Omega]_, opts:OptionsPattern[]] :=
 Module[{\[Lambda], rmin, rmax, BCs, norms, solFuncs, RWRF, m = 0, a=0},
  (* Compute the eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* rmin and rmax *)
  rmin = OptionValue["rmin"];
  rmax = OptionValue["rmax"];
  If[Not[NumericQ[rmin]],
    Message[ReggeWheelerRadial::optx, "rmin" -> rmin];
    Return[$Failed];
  ];
  If[Not[NumericQ[rmax]],
    Message[ReggeWheelerRadial::optx, "rmax" -> rmax];
    Return[$Failed];
  ];
  
  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue[{ReggeWheelerRadial, ReggeWheelerRadialNumericalIntegration}, {opts}, "BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[ReggeWheelerRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Asymptotic normalizations such that we have unit transmission coefficient *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    BCs /. {"In" -> ReggeWheeler`NumericalIntegration`Private`PsiIn[s, l, \[Omega], rmin, rmax]["Psi"],
            "Up" -> ReggeWheeler`NumericalIntegration`Private`PsiUp[s, l, \[Omega], rmin, rmax]["Psi"]};

  (* We only need the normalisations for the specified boundary conditions *)
  norms = BCs /. norms;

  (* Function to construct a ReggeWheelerRadialFunction *)
  RWRF[bc_, ns_, sf_] :=
    ReggeWheelerRadialFunction[s, l, \[Omega],
     Association["s" -> s, "l" -> l, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Method" -> {"NumericalIntegration", "rmin" -> rmin, "rmax" -> rmax},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "RadialFunction" -> sf
     ]
    ];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> RWRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[RWRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*MST Method*)


Options[ReggeWheelerRadialMST] = {
  "RenormalizedAngularMomentum" -> "Monodromy",
  "BoundaryConditions" -> None};


ReggeWheelerRadialMST[s_Integer, l_Integer, \[Omega]_, opts:OptionsPattern[]] :=
 Module[{\[Lambda], \[Nu], BCs, norms, solFuncs, RWRF, m = 0, a=0},
  (* Compute the eigenvalue and renormalized angular momentum *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  \[Nu] = RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method -> OptionValue["RenormalizedAngularMomentum"]];

  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue[{ReggeWheelerRadial, ReggeWheelerRadialMST}, {opts}, "BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[ReggeWheelerRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Compute the asymptotic normalisations *)
  norms = ReggeWheeler`MST`MST`Private`Amplitudes[s, l, m, a, 2\[Omega], \[Nu], \[Lambda]];

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    BCs /. {"In" -> ReggeWheeler`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In"]["Transmission"]],
            "Up" -> ReggeWheeler`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up"]["Transmission"]]};

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms/norms[[All, "Transmission"]];
  norms = BCs /. norms;

  (* Function to construct a ReggeWheelerRadialFunction *)
  RWRF[bc_, ns_, sf_] :=
    ReggeWheelerRadialFunction[s, l, \[Omega],
     Association["s" -> s, "l" -> l, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> \[Nu]},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "RadialFunction" -> sf
     ]
    ];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> RWRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[RWRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*ReggeWheelerRadial*)


SyntaxInformation[ReggeWheelerRadial] =
 {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};


Options[ReggeWheelerRadial] = {
  Method -> "MST",
  "BoundaryConditions" -> {"In", "Up"}
};


ReggeWheelerRadial[s_Integer, l_Integer, \[Omega]_?InexactNumberQ, opts:OptionsPattern[]] :=
 Module[{RWRF, subopts},
  (* All options  except for Method are passed on. Method is a special case where
     only suboptions are passed on, if there are any. *)
  subopts = Cases[{opts}, Except[Method -> _]];
  If[ListQ[OptionValue[Method]],
    subopts = Join[Rest[OptionValue[Method]], subopts];
  ];

  (* Decide which implementation to use *)
  Switch[OptionValue[Method],
    "MST" | {"MST", Rule[_,_]...},
      RWRF = ReggeWheelerRadialMST,
    "NumericalIntegration" | {"NumericalIntegration", ___},
      RWRF = ReggeWheelerRadialNumericalIntegration;,
    _,
      Message[ReggeWheelerRadial::optx, Method -> OptionValue[Method]];
      Return[$Failed];
  ];

  (* Check only supported sub-options have been specified *)  
  If[subopts =!= (subopts = FilterRules[subopts, Options[RWRF]]),
    Message[ReggeWheelerRadial::optx, Method -> OptionValue[Method]];
  ];

  (* Call the chosen implementation *)
  RWRF[s, l, \[Omega], Sequence@@subopts]
];


(* ::Section::Closed:: *)
(*ReggeWheelerRadialFunction*)


Format[ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]] := 
  "ReggeWheelerRadialFunction["<>ToString[s]<>","<>ToString[l]<>","<>ToString[\[Omega]]<>",<<>>]";

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][y_String] /; !MemberQ[{"RadialFunction"}, y] :=
  assoc[y];

ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][r_?NumericQ] :=
  assoc["RadialFunction"][r];

Derivative[n_][ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]][r_?NumericQ] :=
  assoc["RadialFunction"]'[r];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];

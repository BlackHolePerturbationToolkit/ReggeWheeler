(* ::Package:: *)

(* ::Title:: *)
(*ReggeWheeler package*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ReggeWheelerRadial`",
  {
    "ReggeWheeler`NumericalIntegration`",
    "ReggeWheeler`MST`RenormalizedAngularMomentum`",
    "ReggeWheeler`MST`MST`",
    "SpinWeightedSpheroidalHarmonics`"
  }];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{ReggeWheelerRadial, ReggeWheelerRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


ReggeWheelerRadial::usage = "ReggeWheelerRadial[s, l, \[Omega]] computes homogeneous solutions to the Regge Wheeler equation."
ReggeWheelerRadialFunction::usage = "ReggeWheelerRadialFunction[s, l, \[Omega], assoc] is an object representing a homogeneous solution to the Regge Wheeler equation."


(* ::Subsection::Closed:: *)
(*Error Messages*)


ReggeWheelerRadial::precw = "The precision of \[Omega]=`1` is less than WorkingPrecision (`2`).";
ReggeWheelerRadial::optx = "Unknown options in `1`";
ReggeWheelerRadial::dm = "Option `1` is not valid with BoundaryConditions \[RightArrow] `2`.";
ReggeWheelerRadial::sopt = "Option `1` not supported for static (\[Omega]=0) modes.";
ReggeWheelerRadial::hcopt = "Option `1` not supported for HeunC method.";
ReggeWheelerRadialFunction::dmval = "Radius `1` lies outside the computational domain. Results may be incorrect.";
ReggeWheelerRadialFunction::pot = "Invalid potential `1`.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ReggeWheelerRadial*)


(* ::Subsection::Closed:: *)
(*Numerical Integration Method*)


Options[ReggeWheelerRadialNumericalIntegration] = Join[
  {"Domain" -> None},
  FilterRules[Options[NDSolve], Except[WorkingPrecision|AccuracyGoal|PrecisionGoal]]];


domainQ[domain_] := MatchQ[domain, {_?NumericQ, _?NumericQ} | (_?NumericQ) | All];


ReggeWheelerRadialNumericalIntegration[s_Integer, l_Integer, \[Omega]_, BCs_, pot_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], RWRF, norms, ndsolveopts, solFuncs, domains, m = 0, a=0},
  (* Compute the eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* Function to construct a single ReggeWheelerRadialFunction *)
  RWRF[bc_, ns_, sf_, domain_, ndsolveopts___] :=
   Module[{solutionFunction},
    solutionFunction = sf[domain];
    ReggeWheelerRadialFunction[s, l, \[Omega],
     Association["s" -> s, "l" -> l, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Potential" -> pot,
      "Method" -> {"NumericalIntegration", ndsolveopts},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> If[domain === All, {2, \[Infinity]}, First[solutionFunction["Domain"]]],
      "RadialFunction" -> solutionFunction
     ]
    ]
   ];

  (* Domain over which the numerical solution can be evaluated *)
  domains = OptionValue["Domain"];
  If[ListQ[BCs],
    If[!MatchQ[domains, (List|Association)[Rule["In"|"Up",_?domainQ]..]],
      Message[ReggeWheelerRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
    domains = Lookup[domains, BCs, None]; 
    If[!AllTrue[domains, domainQ],
      Message[ReggeWheelerRadial::dm, "Domain" -> OptionValue["Domain"], BCs];
      Return[$Failed];
    ];
  ,
    If[!domainQ[domains],
      Message[ReggeWheelerRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
  ];

  (* Asymptotic normalizations such that we have unit transmission coefficient *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;
  norms = Lookup[norms, BCs];

  (* Solution functions for the specified boundary conditions *)
  ndsolveopts = Sequence@@FilterRules[{opts}, Options[NDSolve]];
  solFuncs =
   <|"In" :> ReggeWheeler`NumericalIntegration`Private`Psi[s, l, \[Omega], "In", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, ndsolveopts],
     "Up" :> ReggeWheeler`NumericalIntegration`Private`Psi[s, l, \[Omega], "Up", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, ndsolveopts]|>;
  solFuncs = Lookup[solFuncs, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> RWRF[#1, #2, #3, #4, ndsolveopts]&, {BCs, norms, solFuncs, domains}]]],
    Return[RWRF[BCs, norms, solFuncs, domains, ndsolveopts]]
  ];
];


(* ::Subsection::Closed:: *)
(*MST Method*)


Options[ReggeWheelerRadialMST] = {
  "RenormalizedAngularMomentum" -> "Monodromy"};


ReggeWheelerRadialMST[s_Integer, l_Integer, \[Omega]_, BCs_, pot_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], \[Nu], norms, solFuncs, RWRF, m = 0, a=0},
  (* Compute the eigenvalue and renormalized angular momentum *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  \[Nu] = RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method -> OptionValue["RenormalizedAngularMomentum"]];

  (* Function to construct a ReggeWheelerRadialFunction *)
  RWRF[bc_, ns_, sf_] :=
    ReggeWheelerRadialFunction[s, l, \[Omega],
     Association["s" -> s, "l" -> l, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Potential" -> pot,
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> \[Nu]},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> {2, \[Infinity]}, "RadialFunction" -> sf
     ]
    ];

  (* Compute the asymptotic normalisations *)
  norms = ReggeWheeler`MST`MST`Private`Amplitudes[s, l, m, a, 2\[Omega], \[Nu], \[Lambda], {wp, prec, acc}];

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    <|"In" :> ReggeWheeler`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In"]["Transmission"], {wp, prec, acc}],
      "Up" :> ReggeWheeler`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up"]["Transmission"], {wp, prec, acc}]|>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms[[All, {"Transmission"}]]/norms[[All, "Transmission"]];
  norms = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> RWRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[RWRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*HeunC*)


Options[ReggeWheelerRadialHeunC] = {};


ReggeWheelerRadialHeunC[s_Integer, l_Integer, \[Omega]_, BCs_, pot_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], \[Nu], norms, solFuncs, RWRF, m = 0, a=0},
  (* Compute the eigenvalue and renormalized angular momentum *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* Function to construct a ReggeWheelerRadialFunction *)
  RWRF[bc_, ns_, sf_] :=
    ReggeWheelerRadialFunction[s, l, \[Omega],
     Association["s" -> s, "l" -> l, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Potential" -> pot,
      "Method" -> {"HeunC"},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> {2, \[Infinity]}, "RadialFunction" -> sf
     ]
    ];

  (* Compute the asymptotic normalisations *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    <|"In" :> Function[{r}, (r/2)^(s + 1) E^(-I \[Omega] (r + 2 Log[r/2 - 1])) HeunC[l + l^2 - (1 + s) (s - 4 I \[Omega]), 4 I (1 + s) \[Omega], 1 - 4 I \[Omega], 1 + 2 s, 4 I \[Omega], 1 - r/2]],
      "Up" :> Function[{r}, $Failed]|>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms[[All, {"Transmission"}]]/norms[[All, "Transmission"]];
  norms = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> RWRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[RWRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*Static modes*)


ReggeWheelerRadialStatic[s_Integer, l_Integer, \[Omega]_, BCs_, pot_] :=
 Module[{\[Lambda], norms, solFuncs, RWRF, m = 0, a=0},
  (* Compute the eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, 0];

  (* Function to construct a ReggeWheelerRadialFunction *)
  RWRF[bc_, ns_, sf_] :=
    ReggeWheelerRadialFunction[s, l, \[Omega],
     Association["s" -> s, "l" -> l, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Potential" -> pot,
      "Method" -> {"Static"},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> {2, \[Infinity]}, "RadialFunction" -> sf
     ]
    ];

  (* Compute the asymptotic normalisations *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    <|"In" :> Function[{r},(r/2)^(-l)*Hypergeometric2F1[l+s+1,l-s+1,1,(r-2)/r]],
      "Up" :> Function[{r},(r/2)^(-l)*Hypergeometric2F1[l+s+1,l-s+1,2*(l+1),2/r]]|>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms[[All, {"Transmission"}]]/norms[[All, "Transmission"]];
  norms = Lookup[norms, BCs];

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
  Method -> Automatic,
  "BoundaryConditions" -> {"In", "Up"},
  "Potential" -> "ReggeWheeler",
  WorkingPrecision -> Automatic,
  PrecisionGoal -> Automatic,
  AccuracyGoal -> Automatic
};


(* ::Subsubsection::Closed:: *)
(*Static modes*)


ReggeWheelerRadial[s_Integer, l_Integer, \[Omega]_, opts:OptionsPattern[]] /; \[Omega] == 0 :=
 Module[{pot, BCs, wp, prec, acc},
  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue["BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[ReggeWheelerRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Potential *)
  pot = OptionValue["Potential"];

  (* Options are not supported for static modes *)
  Do[
    If[OptionValue[opt] =!= Automatic, Message[ReggeWheelerRadial::sopt, opt]];,
    {opt, {Method, WorkingPrecision, PrecisionGoal, AccuracyGoal}}
  ];

  (* Call the chosen implementation *)
  ReggeWheelerRadialStatic[s, l, \[Omega], BCs, pot]
]


(* ::Subsubsection::Closed:: *)
(*Non-static modes*)


ReggeWheelerRadial[s_Integer, l_Integer, \[Omega]_?InexactNumberQ, opts:OptionsPattern[]] :=
 Module[{RWRF, subopts, pot, BCs, wp, prec, acc},
  (* Extract suboptions from Method to be passed on. *)
  If[ListQ[OptionValue[Method]],
    subopts = Rest[OptionValue[Method]];,
    subopts = {};
  ];

  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue["BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[ReggeWheelerRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Potential *)
  pot = OptionValue["Potential"];

  (* Options associated with precision and accuracy *)
  {wp, prec, acc} = OptionValue[{WorkingPrecision, PrecisionGoal, AccuracyGoal}];
  If[wp === Automatic, wp = Precision[\[Omega]]];
  If[prec === Automatic, prec = wp / 2];
  If[acc === Automatic, acc = wp / 2];
  If[Precision[\[Omega]] < wp, Message[ReggeWheelerRadial::precw, \[Omega], wp]];

  (* Decide which implementation to use *)
  Switch[OptionValue[Method],
    Automatic,
      RWRF = ReggeWheelerRadialMST,
    "MST" | {"MST", OptionsPattern[ReggeWheelerRadialMST]},
      RWRF = ReggeWheelerRadialMST,
    "NumericalIntegration" | {"NumericalIntegration", OptionsPattern[ReggeWheelerRadialNumericalIntegration]},
      RWRF = ReggeWheelerRadialNumericalIntegration;,
    "HeunC",
      (* Some options are not supported for the HeunC method modes *)
      Do[
        If[OptionValue[opt] =!= Automatic, Message[ReggeWheelerRadial::hcopt, opt]];,
        {opt, {WorkingPrecision, PrecisionGoal, AccuracyGoal}}
      ];
      RWRF = ReggeWheelerRadialHeunC;,
    _,
      Message[ReggeWheelerRadial::optx, Method -> OptionValue[Method]];
      Return[$Failed];
  ];

  (* Check only supported sub-options have been specified *)  
  If[subopts =!= (subopts = FilterRules[subopts, Options[RWRF]]),
    Message[ReggeWheelerRadial::optx, Method -> OptionValue[Method]];
  ];

  (* Call the chosen implementation *)
  RWRF[s, l, \[Omega], BCs, pot, {wp, prec, acc}, Sequence@@subopts]
];


(* ::Section::Closed:: *)
(*ReggeWheelerRadialFunction*)


(* ::Subsection::Closed:: *)
(*Output format*)


(* ::Subsubsection::Closed:: *)
(*Icons*)


icons = <|
 "In" -> Graphics[{
         Line[{{0,1/2},{1/2,1},{1,1/2},{1/2,0},{0,1/2}}],
         Line[{{3/4,1/4},{1/2,1/2}}],
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{1/4,3/4}}]]},
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{3/4,3/4}}]]}},
         Background -> White,
         ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]],
 "Up" -> Graphics[{
         Line[{{0,1/2},{1/2,1},{1,1/2},{1/2,0},{0,1/2}}],
         Line[{{1/4,1/4},{1/2,1/2}}],
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{1/4,3/4}}]]},
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{3/4,3/4}}]]}},
         Background -> White,
         ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]]
|>;


(* ::Subsubsection::Closed:: *)
(*Formatting of ReggeWheelerRadialFunction*)


ReggeWheelerRadialFunction /:
 MakeBoxes[rwrf:ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", s}], "  ",
                  BoxForm`SummaryItem[{"l: ", l}], "  ",
                  BoxForm`SummaryItem[{"\[Omega]: ", \[Omega]}]}],
             BoxForm`SummaryItem[{"Domain: ", assoc["Domain"]}],
             BoxForm`SummaryItem[{"Boundary Conditions: " , assoc["BoundaryConditions"]}]};
  extended = {BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Transmission Amplitude: ", assoc["Amplitudes", "Transmission"]}],
              BoxForm`SummaryItem[{"Incidence Amplitude: ", Lookup[assoc["Amplitudes"], "Incidence", Missing]}],
              BoxForm`SummaryItem[{"Reflection Amplitude: ", Lookup[assoc["Amplitudes"], "Reflection", Missing]}],
              BoxForm`SummaryItem[{"Method: ", First[assoc["Method"]]}],
              BoxForm`SummaryItem[{"Method options: ",Column[Rest[assoc["Method"]]]}]};
  BoxForm`ArrangeSummaryBox[
    ReggeWheelerRadialFunction,
    rwrf,
    Lookup[icons, assoc["BoundaryConditions"], None],
    summary,
    extended,
    form]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][y_String] /; !MemberQ[{"RadialFunction"}, y] :=
  assoc[y];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


SetAttributes[ReggeWheelerRadialFunction, {NumericFunction}];


outsideDomainQ[r_, rmin_, rmax_] := Min[r]<rmin || Max[r]>rmax;


ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax, \[Lambda], sign, R},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[ReggeWheelerRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
  ];
  R = assoc["RadialFunction"];
  Quiet[Switch[assoc["Potential"],
    "ReggeWheeler",
    R[r],
    "Zerilli",
    sign = Switch[assoc["BoundaryConditions"], "In", -1, "Out", 1, _, Indeterminate];
    \[Lambda] = assoc["Eigenvalue"];
    1/(\[Lambda]^2 + \[Lambda] + sign 3 I \[Omega]) ((\[Lambda]^2+\[Lambda]+(9(r-2))/(r^2 (3+\[Lambda] r)))R[r]+3(1-2/r)R'[r]),
    _,
    Message[ReggeWheelerRadialFunction::pot, assoc["Potential"]]
    ]
  , InterpolatingFunction::dmval]
 ];


Derivative[n_][ReggeWheelerRadialFunction[s_, l_, \[Omega]_, assoc_]][r0:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax, \[Lambda], sign, R, r},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r0, rmin, rmax],
    Message[ReggeWheelerRadialFunction::dmval, #]& /@ Select[Flatten[{r0}], outsideDomainQ[#, rmin, rmax]&];
  ];
  Quiet[Switch[assoc["Potential"],
    "ReggeWheeler",
    Derivative[n][assoc["RadialFunction"]][r0],
    "Zerilli",
    sign = Switch[assoc["BoundaryConditions"], "In", -1, "Out", 1, _, Indeterminate];
    \[Lambda] = assoc["Eigenvalue"];
    (* FIXME: we could reduce this using the Zerilli equation *)
    Collect[D[1/(\[Lambda]^2 + \[Lambda] + sign 3 I \[Omega]) ((\[Lambda]^2+\[Lambda]+(9(r-2))/(r^2 (3+\[Lambda] r)))R[r]+3(1-2/r)R'[r]),{r,n}], {Derivative[_][R][r], R[r]}] /. R -> assoc["RadialFunction"] /. r -> r0,
    _,
    Message[ReggeWheelerRadialFunction::pot, assoc["Potential"]]
    ]
  , InterpolatingFunction::dmval]
 ];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{ReggeWheelerRadial, ReggeWheelerRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[]
EndPackage[];

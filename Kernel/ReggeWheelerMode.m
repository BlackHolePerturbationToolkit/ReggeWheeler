(* ::Package:: *)

(* ::Title:: *)
(*ReggeWheelerMode*)


(* ::Section:: *)
(*Create Package*)


(* ::Subsection:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ReggeWheelerMode`",
  {"ReggeWheeler`ReggeWheelerSource`",
   "ReggeWheeler`ReggeWheelerRadial`",
   "ReggeWheeler`ConvolveSource`",
   "KerrGeodesics`KerrGeoOrbit`",
   "KerrGeodesics`OrbitalFrequencies`",
   "SpinWeightedSpheroidalHarmonics`"}
];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{ReggeWheelerMode, ReggeWheelerPointParticleMode}, {Protected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


ReggeWheelerMode::usage = "ReggeWheelerMode[assoc] is an object which represents a Regge Wheeler mode.";
ReggeWheelerPointParticleMode::usage = "ReggeWheelerPointParticleMode[s, l, m, n, orbit] produces a "<>
 "ReggeWheelerMode representing a solution to the Regge-Wheeler equation with a point particle source.";


(* ::Subsection::Closed:: *)
(*Error Messages*)


ReggeWheelerPointParticleMode::nospin = "Regge-Wheeler perturbations are only for Schwarzschild black holes but spin `1` is not zero."


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*ReggeWheelerPointParticleMode*)


SyntaxInformation[ReggeWheelerPointParticleMode] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[ReggeWheelerPointParticleMode] = {Method->Automatic};


ReggeWheelerPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]] /; AllTrue[orbit["Frequencies"], InexactNumberQ] :=
 Module[{source, assoc, radialopts, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z,hypopts},
  If[orbit["a"] != 0,
    Message[ReggeWheelerPointParticleMode::nospin, orbit["a"]];
    Return[$Failed];
  ];
  
  If[OptionValue["Method"] == "Hyperboloidal", Return[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal[s,l,m,n,orbit,OptionValue["Method"][[2;;]]]]];
  
  Switch[OptionValue["Method"],
       "Hyperboloidal", Return[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal[s,l,m,n,orbit]],
       {"Hyperboloidal",OptionsPattern[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal]}, 
       hypopts = Sequence@@FilterRules[{OptionValue["Method"][[2;;]]}, Options[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal]];
       Return[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal[s,l,m,n,orbit,hypopts]]
       ];
       

  (*{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];*) (*This gives Mino frequencies, need BL frequencies*)
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = Values[KerrGeoFrequencies[orbit["a"], orbit["p"], orbit["e"], orbit["Inclination"]]];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r;

  source = ReggeWheeler`ReggeWheelerSource`Private`ReggeWheelerPointParticleSource[s, l, m, orbit];

  radialopts = Sequence@@FilterRules[{opts}, Options[ReggeWheelerRadial]]; (*May need to remove the "Potential" option if it's in this list!*)
  If[EvenQ[l+m], (*Switch for parity of the homogeneous solution*)
      R = ReggeWheelerRadial[s, l, \[Omega], "Potential"->"Zerilli", radialopts];
  ,
      R = ReggeWheelerRadial[s, l, \[Omega], "Potential"->"ReggeWheeler", radialopts];
  ];
  S = SpinWeightedSpheroidalHarmonicS[s, l, m, 0];
  Z = ReggeWheeler`ConvolveSource`Private`ConvolveSource[R, S, source];

  assoc = <| "s" -> s, 
             "l" -> l,
             "m" -> m,
             "\[Omega]" -> \[Omega],
             "Eigenvalue" -> R["In"]["Eigenvalue"],
             "Type" -> {"PointParticleCircular", "Radius" -> orbit["p"]},
             "RadialFunctions" -> R,
             "AngularFunction" -> S,
             "Amplitudes" -> Z,
             "Method" -> R["In"]["Method"]
           |>;

  ReggeWheelerMode[assoc]
]


(* ::Section:: *)
(*ReggeWheelerMode*)


(* ::Subsection::Closed:: *)
(*Output format*)


ReggeWheelerMode /:
 MakeBoxes[rwm:ReggeWheelerMode[assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", assoc["s"]}], "  ",
                  BoxForm`SummaryItem[{"l: ", assoc["l"]}], "  ",
                  BoxForm`SummaryItem[{"m: ", assoc["m"]}], "  ",
                  BoxForm`SummaryItem[{"\[Omega]: ", assoc["\[Omega]"]}]}],
             BoxForm`SummaryItem[{"Type: ", First[assoc["Type"]]}]};
  extended = {BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Amplitude at \[ScriptCapitalI]: ", assoc["Amplitudes"]["\[ScriptCapitalI]"]}],
              BoxForm`SummaryItem[{"Amplitude at \[ScriptCapitalH]: ", assoc["Amplitudes"]["\[ScriptCapitalH]"]}],
              BoxForm`SummaryItem[{"Type details: ", Column[Rest[assoc["Type"]]]}]};
  BoxForm`ArrangeSummaryBox[
    ReggeWheelerMode,
    rwm,
    None,
    summary,
    extended,
    form
  ]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


ReggeWheelerMode[assoc_]["EnergyFlux"] := EnergyFlux[ReggeWheelerMode[assoc]];


ReggeWheelerMode[assoc_]["Fluxes"] := <|"Energy" -> ReggeWheelerMode[assoc]["EnergyFlux"], "AngularMomentum" -> ReggeWheelerMode[assoc]["AngularMomentumFlux"]|>;


ReggeWheelerMode[assoc_]["AngularMomentumFlux"] := AngularMomentumFlux[ReggeWheelerMode[assoc]];


ReggeWheelerMode[assoc_][string_] := assoc[string];


Keys[m_ReggeWheelerMode] ^:= Keys[m[[1]]]


(* ::Section::Closed:: *)
(*Fluxes*)


(* ::Subsection::Closed:: *)
(*Energy Flux*)


EnergyFlux[mode_ReggeWheelerMode] :=
 Module[{l, m, \[Omega], Z, FluxInf, FluxH,\[Xi]},
  l = mode["l"];
  m = mode["m"];
  \[Omega] = mode["\[Omega]"];
  Z = mode["Amplitudes"];
  
  If[mode["Method"][[1]] == "Hyperboloidal",
      \[Xi] = -I \[Omega] 4;
  	FluxInf = ((l+2)!/(l-2)!)1/(256\[Pi]*16)Abs[\[Xi] Z[[1]]]^2;
	  FluxH = ((l+2)!/(l-2)!)1/(256\[Pi]*16)Abs[\[Xi] Z[[2]]]^2;
      ,
      FluxInf = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1))*Abs[\[Omega]*Z["\[ScriptCapitalI]"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2))*Abs[\[Omega]*Z["\[ScriptCapitalI]"]]^2/(16*Pi)];
      FluxH   = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1))*Abs[\[Omega]*Z["\[ScriptCapitalH]"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2))*Abs[\[Omega]*Z["\[ScriptCapitalH]"]]^2/(16*Pi)];
   ];
  
  <| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxH |>
];


(* ::Subsection::Closed:: *)
(*Angular Momentum Flux*)


AngularMomentumFlux[mode_ReggeWheelerMode] :=
 Module[{l, m, \[Omega], Z, FluxInf, FluxH,\[Xi]},
  l = mode["l"];
  m = mode["m"];
  \[Omega] = mode["\[Omega]"];
  Z = mode["Amplitudes"];

 If[mode["Method"][[1]] == "Hyperboloidal",
 	\[Xi] = -I \[Omega] 4;
	 FluxInf = I*m*\[Xi]((l+2)!/(l-2)!)1/(64\[Pi]*16)Abs[Z[[1]]]^2;
	 FluxH = I*m*\[Xi]((l+2)!/(l-2)!)1/(64\[Pi]*16)Abs[Z[[2]]]^2;
	 ,
     FluxInf = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1)) m \[Omega] Abs[Z["\[ScriptCapitalI]"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2)) m \[Omega] Abs[Z["\[ScriptCapitalI]"]]^2/(16*Pi)];
     FluxH   = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1)) m \[Omega] Abs[Z["\[ScriptCapitalH]"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2)) m \[Omega] Abs[Z["\[ScriptCapitalH]"]]^2/(16*Pi)];
  ];
  
  <| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxH |>
];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{ReggeWheelerMode, ReggeWheelerPointParticleMode}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];

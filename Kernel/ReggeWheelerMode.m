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


(* ::Subsection:: *)
(*Error Messages*)


ReggeWheelerPointParticleMode::nospin = "Regge-Wheeler perturbations are only for Schwarzschild black holes but spin `1` is not zero."
ReggeWheelerPointParticleMode::eccentricity = "The hyperboloidal package does not currently accept eccentric orbits. Please set eccentricity ('e') to zero."
ReggeWheelerPointParticleMode::spin2field = "The hyperboloidal package currently only works for spin = 2 fields, but the fluxes and radial fns. are correct for spin = -2. Please set field spin ('s') to two."
ReggeWheelerPointParticleMode::inclination = "The hyperboloidal package currently only works for orbits in the equatorial plane. Please set orbital inclination ('x') to one."
ReggeWheelerPointParticleMode::eccentricitymode = "The hyperboloidal package currently only works for circular orbit modes ('m'). Please set the eccentricity mode ('n') to zero."


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
  
  
  Switch[OptionValue["Method"],
       "Hyperboloidal",	   
       If[orbit["e"] != 0,
		    Message[ReggeWheelerPointParticleMode::eccentricity];
		    Return[$Failed];
	   ];
			
	   If[s != 2,
		    Message[ReggeWheelerPointParticleMode::spin2field];
		    Return[$Failed];
	   ];
			
	   If[orbit["Inclination"] != 1,
		    Message[ReggeWheelerPointParticleMode::inclination];
		    Return[$Failed];
	   ];
			
	   If[n != 0,
		    Message[ReggeWheelerPointParticleMode::eccentricitymode];
		    Return[$Failed];
	   ];
	    Return[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal[s,l,m,n,orbit]],
       
       {"Hyperboloidal",OptionsPattern[ReggeWheeler`Hyperboloidal`Private`ReggeWheelerHyperboloidal]},
       If[orbit["e"] != 0,
		    Message[ReggeWheelerPointParticleMode::eccentricity];
		    Return[$Failed];
	   ];
			
	   If[s != 2,
		    Message[ReggeWheelerPointParticleMode::spin2field];
		    Return[$Failed];
	   ];
			
	   If[orbit["Inclination"] != 1,
		    Message[ReggeWheelerPointParticleMode::inclination];
		    Return[$Failed];
	   ];
			
	   If[n != 0,
		    Message[ReggeWheelerPointParticleMode::eccentricitymode];
		    Return[$Failed];
	   ];
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


(* ::Section:: *)
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
  	FluxInf = ((l+2)!/(l-2)!)1/(256\[Pi]*16)Abs[-I \[Omega] 4 Z[[1]]]^2;
	  FluxH = ((l+2)!/(l-2)!)1/(256\[Pi]*16)Abs[-I \[Omega] 4 Z[[2]]]^2;
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
	 FluxInf = I*m*(-I \[Omega] 4)((l+2)!/(l-2)!)1/(64\[Pi]*16)Abs[Z[[1]]]^2;
	 FluxH = I*m*(-I \[Omega] 4)((l+2)!/(l-2)!)1/(64\[Pi]*16)Abs[Z[[2]]]^2;
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

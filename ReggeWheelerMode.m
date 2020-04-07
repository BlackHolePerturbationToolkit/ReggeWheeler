(* ::Package:: *)

(* ::Title:: *)
(*ReggeWheelerMode*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ReggeWheelerMode`",
  {"ReggeWheeler`ReggeWheelerSource`",
   "ReggeWheeler`ReggeWheelerRadial`",
   "ReggeWheeler`ConvolveSource`",
   "KerrGeodesics`KerrGeoOrbit`",
   "KerrGeodesics`OrbitalFrequencies`"}
];


(* ::Subsection::Closed:: *)
(*Usage messages*)


ReggeWheelerPointParticleMode::usage = "ReggeWheelerPointParticleMode[s, l, m, n, orbit] solves the Regge Wheeler equation with a point particle source.";
ReggeWheelerMode::usage = "ReggeWheelerMode[assoc] is an object which represents a Regge Wheeler mode.";
EnergyFlux::usage = "EnergyFlux[mode] computes the flux of energy radiated in the given mode.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ReggeWheelerPointParticleMode*)


(*defined currently for circular orbits*)
ReggeWheelerPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction] :=
	Module[{source},
	
		If[!PossibleZeroQ[orbit["a"]], Print["Regge-Wheeler perturbations are only for Schwarzschild (a=0) black holes"]; Return[];];
	
		source = ReggeWheelerPointParticleSource[s, l, m, orbit];
		Return[ReggeWheelerPointParticleMode[s, l, m, n, orbit, source]];
	];


ReggeWheelerPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, source_ReggeWheelerSourceObject] :=
 Module[{assoc, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z},
  (*{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];*) (*This gives Mino frequencies, need BL frequencies*)
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = Values[KerrGeoFrequencies[orbit["a"], orbit["p"], orbit["e"], orbit["Inclination"]]];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r;

  R = ReggeWheelerRadial[s, l, \[Omega]];
  

  Z = ReggeWheeler`ConvolveSource`Private`ConvolveSource[R, source];

  assoc = <| "s" -> s, 
             "l" -> l,
             "m" -> m,
             "n" -> n,
             "Type" -> "PointParticleCircular",
             "Homogeneous" -> R,
             "Particular" -> Z
           |>;

  ReggeWheelerMode[assoc]
]


(* ::Section::Closed:: *)
(*ReggeWheelerMode*)


(* ::Subsection::Closed:: *)
(*Output format*)


Format[ReggeWheelerModeObject[assoc_]] := "ReggeWheelerModeObject["<>ToString[assoc["s"]]<>","<>ToString[assoc["l"]]<>","<>ToString[assoc["m"]]<>","<>ToString[assoc["n"]]<>","<>"<<>>]";


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


ReggeWheelerMode[assoc_][string_] := assoc[string];


(* ::Section::Closed:: *)
(*Fluxes*)


(* ::Subsection::Closed:: *)
(*Energy Flux*)


EnergyFlux[mode_ReggeWheelerMode] :=
 Module[{l, m, \[Omega], Z, FluxInf, FluxH},
  l = mode["l"];
  m = mode["m"];
  \[Omega] = mode["\[Omega]"];
  Z = mode["Amplitudes"];

  FluxInf = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1))*Abs[\[Omega]*Z["\[ScriptCapitalI]"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2))*Abs[\[Omega]*Z["\[ScriptCapitalI]"]]^2/(16*Pi)];
  FluxH   = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1))*Abs[\[Omega]*Z["\[ScriptCapitalH]"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2))*Abs[\[Omega]*Z["\[ScriptCapitalH]"]]^2/(16*Pi)];
  
  <| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxH |>
];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];

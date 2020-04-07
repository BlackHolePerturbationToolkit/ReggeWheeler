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


ReggeWheelerModeObject::usage = "ReggeWheelerModeObject[assoc] an object which contains a Regge Wheeler mode.";
ReggeWheelerPointParticleMode::usage = "ReggeWheelerPointParticleMode[s, l, m, n, orbit] solves the Regge Wheeler equation with a point particle source.";


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
 Module[{assoc, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z, FluxH, FluxInf, Fluxes},
  (*{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];*) (*This gives Mino frequencies, need BL frequencies*)
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = Values[KerrGeoFrequencies[orbit["a"], orbit["p"], orbit["e"], orbit["Inclination"]]];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r;

  R = ReggeWheelerRadial[s, l, \[Omega]];
  

  Z = ReggeWheeler`ConvolveSource`Private`ConvolveSource[R, source];
  
  (*currently fluxes for circular orbit*)
   FluxInf = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZInf"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZInf"]]^2/(16*Pi)];
   FluxH   = If[EvenQ[l+m], (l-1)*(l+2)/(l*(l+1))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZHor"]]^2/(4*Pi), (l*(l+1))/((l-1)*(l+2))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZHor"]]^2/(16*Pi)];
  
  
  Fluxes = <| "FluxInf" -> FluxInf, "FluxHor" -> FluxH, "FluxTotal" -> FluxH + FluxInf |>;

  assoc = <| "s" -> s, 
             "l" -> l,
             "m" -> m,
             "n" -> n,
             "Type" -> "PointParticleCircular",
             "Homogeneous" -> R,
             "Particular" -> Z,
             "Fluxes" -> Fluxes
           |>;

  Return[ReggeWheelerModeObject[assoc]];
]


(* ::Section::Closed:: *)
(*ReggeWheelerMode*)


(* ::Subsection::Closed:: *)
(*Output format*)


Format[ReggeWheelerModeObject[assoc_]] := "ReggeWheelerModeObject["<>ToString[assoc["s"]]<>","<>ToString[assoc["l"]]<>","<>ToString[assoc["m"]]<>","<>ToString[assoc["n"]]<>","<>"<<>>]";


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


ReggeWheelerModeObject[assoc_][string_] := assoc[string]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];

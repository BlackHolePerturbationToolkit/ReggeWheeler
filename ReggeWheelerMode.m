(* ::Package:: *)

BeginPackage["ReggeWheeler`ReggeWheelerMode`",
  {"ReggeWheeler`ReggeWheelerSource`",
   "ReggeWheeler`ReggeWheelerRadial`",
   "ReggeWheeler`ConvolveSource`",
   "KerrGeodesics`KerrGeoOrbit`"}
];

ReggeWheelerModeObject::usage = "ReggeWheelerModeObject[assoc] an object which contains a Regge Wheeler mode.";

ReggeWheelerPointParticleMode::usage = "ReggeWheelerPointParticleMode[s, l, m, n, k, orbit] solves the Regge Wheeler equation with a point particle source.";

Begin["`Private`"];


(*defined currently for circular orbits*)
ReggeWheelerPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction] :=
	Module[{source},
		source = ReggeWheelerPointParticleSource[s, l, m, orbit];
		Return[ReggeWheelerPointParticleMode[s, l, m, n, k, orbit, source]];
	];


ReggeWheelerPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction, source_ReggeWheelerSourceObject] :=
 Module[{assoc, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z, Fluxes},
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r;

  R = ReggeWheelerRadial[s, l, \[Omega]];

  Z = ReggeWheeler`ConvolveSource`Private`ConvolveSource[R, source];
  
  (*currently fluxes for circular orbit*)
	Fluxes = If[EvenQ[l+m],
	<|
		"FluxInf" -> (l-1)*(l+2)/(l*(l+1))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZInf"]]^2/(4*Pi),
		"FluxHor" -> (l-1)*(l+2)/(l*(l+1))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZHor"]]^2/(4*Pi)
	|>,
	<|
		"FluxInf" -> (l*(l+1))/((l-1)*(l+2))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZInf"]]^2/(16*Pi),
		"FluxHor" -> (l*(l+1))/((l-1)*(l+2))*Abs[m*\[CapitalOmega]\[Phi]*Z["ZHor"]]^2/(16*Pi)
	|>
	];

  assoc = <| "l" -> l,
             "m" -> m,
             "n" -> n,
             "k" -> k,
             "Type" -> "PointParticleCircular",
             "Homogeneous" -> R,
             "Particular" -> Z,
             "Fluxes" -> Fluxes
           |>;

  Return[ReggeWheelerModeObject[assoc]];
]


Format[ReggeWheelerModeObject[assoc_]] := "ReggeWheelerModeObject[<<>>]";

ReggeWheelerModeObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];

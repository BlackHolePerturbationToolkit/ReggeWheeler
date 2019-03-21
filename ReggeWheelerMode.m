BeginPackage["ReggeWheeler`ReggeWheelerMode`",
  {"ReggeWheeler`ReggeWheelerSource`",
   "ReggeWheeler`ReggeWheelerRadial`",
   "ReggeWheeler`ConvolveSource`"}
];

ReggeWheelerModeObject::usage = "ReggeWheelerModeObject[assoc] an object which contains a Regge Wheeler mode.";

ReggeWheelerPointParticleMode::usage = "ReggeWheelerPointParticleMode[s, l, m, n, orbit] solves the Regge Wheeler equation with a point particle source.";

Begin["`Private`"];


ReggeWheelerPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, source_ReggeWheelerSourceObject] :=
 Module[{assoc, R, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z},
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r;

  R = ReggeWheelerRadial[s, l, \[Omega]];

  Z = ReggeWheeler`ConvolveSource`Private`ConvolveSource[R, source];

  assoc = <| "l" -> l,
             "m" -> m,
             "n" -> n,
             "k" -> k,
             "Type" -> "PointParticle",
             "Radial" -> R,
             "Amplitudes" -> Z
           |>;

  ReggeWheelerModeObject[assoc]
]


Format[ReggeWheelerModeObject[assoc_]] := "ReggeWheelerModeObject[<<>>]";

ReggeWheelerModeObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];

BeginPackage["ReggeWheeler`ReggeWheelerSource`"];

ReggeWheelerSourceObject::usage = "ReggeWheelerSourceObject[assoc] an object which contains a Regge Wheeler source."

ReggeWheelerPointParticleSource::usage = "ReggeWheelerPointParticleSource[s, orbit] Point particle source for the Regge Wheeler equation."

Begin["`Private`"];


ReggeWheelerPointParticleSource[variable_, orbit_] :=
 Module[{assoc, source},
  assoc = <| "Variable" -> variable, "sourceType" -> "PointParticle" |>;
  ReggeWheelerSourceObject[assoc]
];

Format[ReggeWheelerSourceObject[assoc_]] := "ReggeWheelerSourceObject[<<>>]";

ReggeWheelerSourceObject[assoc_][string_] := assoc[string];


End[];
EndPackage[];

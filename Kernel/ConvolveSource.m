(* ::Package:: *)

(* ::Title:: *)
(*ConvolveSource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ConvolveSource`"];


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ConvolveSource*)


ConvolveSource[RF_, SF_, SO_] :=Module[{s},
		s = RF["In"]["s"];
		
		If[Abs[s] == 2, 
			If[SO["type"]=="PointParticleCircular", Return[ConvolvePointParticleSourceCircular[s,RF,SO]], Print["Only circular orbit sources are implemented"]];
			,Print["Only the |s|=2 source is currently implemented"];
		];
	]


ConvolvePointParticleSourceCircular[(-2|2),RF_,SO_]:=
	Module[{r0,PsiIn,dPsiIn,PsiUp,dPsiUp,Wronskian,deltadPsidr,deltaPsi,ZIn,ZUp},

		r0=SO["r0"];

		PsiIn=RF["In"][r0];
		dPsiIn=RF["In"]'[r0];
		PsiUp=RF["Up"][r0];
		dPsiUp=RF["Up"]'[r0];
		
		Wronskian = PsiIn*dPsiUp - PsiUp*dPsiIn;
		deltaPsi=SO["deltaPsi"];
		deltadPsidr=SO["deltadPsidr"];
		ZIn = (PsiUp*deltadPsidr - deltaPsi*dPsiUp)/Wronskian;
		ZUp = (PsiIn*deltadPsidr - deltaPsi*dPsiIn)/Wronskian;
		<|"\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn|>
];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];

(* ::Package:: *)

(* ::Title:: *)
(*ConvolveSource*)


(* ::Section:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ConvolveSource`",{"KerrGeodesics`KerrGeoOrbit`","KerrGeodesics`OrbitalFrequencies`","KerrGeodesics`FourVelocity`"}];


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*ConvolveSource*)


(* ::Subsection::Closed:: *)
(*|s|=2 point particle on a circular orbit*)


ConvolveSource[RF_, SF_, SO_] :=Module[{s},
		s = RF["In"]["s"];
		
		If[Abs[s] == 2, 
			Which[SO["type"]=="PointParticleCircular", Return[ConvolvePointParticleSourceCircular[s,RF,SO]],SO["type"]=="PointParticleEccentric",Return[ConvolvePointParticleSourceEccentric[s,RF,SO]],True, Print["Only circular and eccentric equatorial orbit sources are implemented"]];
			,Print["Only the |s|=2 source is currently implemented"];
		];
	];


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


(* ::Subsection:: *)
(*|s|=2 point particle on a eccentric orbit*)


ConvolvePointParticleSourceEccentric[(-2|2),RF_,SO_]:=
	Module[{M=1,rp,tp,\[Phi]p,dtd\[Chi],\[Omega],\[CapitalOmega]r,p,e,x,PsiIn,dPsiIn,PsiUp,dPsiUp,W,ZInToint,ZUpToint,ZIn,ZUp, OrbitDarwin, FourVelDarwin,ut,ur,u\[Phi],l,m,n,Y,F,G,\[Lambda]},
	
		(*Orbital parameters*)
		\[Omega]=SO["\[Omega]"];
		p=SO["p"];
		e=SO["e"];
		\[CapitalOmega]r=SO["\[CapitalOmega]r"];
		x=SO["x"];
		
		(*Orbit coordinates from KerrGeodesics*)
		OrbitDarwin=KerrGeoOrbit[0,p,e,x,Parametrization->"Darwin"]["Trajectory"];
		rp=OrbitDarwin[[2]];
		tp=tp=OrbitDarwin[[1]];
		\[Phi]p=OrbitDarwin[[4]];
		
		(*Orbit four velocity (upper indices) from KerrGeodesics*)
		FourVelDarwin=KerrGeoFourVelocity[0,p,e,x,Parametrization->"Darwin"];
		ut=FourVelDarwin["\!\(\*SuperscriptBox[\(u\), \(t\)]\)"];
		ur=FourVelDarwin["\!\(\*SuperscriptBox[\(u\), \(r\)]\)"];
		u\[Phi]=FourVelDarwin["\!\(\*SuperscriptBox[\(u\), \(\[Phi]\)]\)"];
		
		(*Homogeneous Solutions*)
		PsiIn = Function[\[Chi],Evaluate[RF["In"][rp[\[Chi]]]]];
		dPsiIn = Function[\[Chi],Evaluate[RF["In"]'[rp[\[Chi]]]]];
		PsiUp = Function[\[Chi],Evaluate[RF["Up"][rp[\[Chi]]]]];
		dPsiUp = Function[\[Chi],Evaluate[RF["Up"]'[rp[\[Chi]]]]];
		
		(*Wronskian (in rstar) is a constant, evaluate at any radius e.g. rp[\[Chi]=0].*)
		W = (1-2 (2M)/rp[0])(PsiIn[0]*dPsiUp[0] - PsiUp[0]*dPsiIn[0]);
		
		(*Check all expressions below, check W should be rstar Wronskian*)
		(*RWZ Source*)
		l=SO["l"];
		m=SO["m"];
		n=SO["n"];
		\[Lambda]=((l+2)(l-1))/2;
		
		(*Principle harmonic function in each sector*)
		Y=If[OddQ[l+m],SphericalHarmonicY[l,m+1,\[Pi]/2,0],SphericalHarmonicY[l,m,\[Pi]/2,0]];
		
		(*Functions in the Jump conditions a-la Hopper+Evans*)
		F=If[OddQ[l+m],Function[\[Chi],Evaluate[(16 E^(-I m \[Phi]p[\[Chi]]) Sqrt[(l-m) (1+l+m)] \[Pi] (-((rp[\[Chi]]^2 ur[\[Chi]]^2)/ut[\[Chi]])+(-2 M+rp[\[Chi]])^2 ut[\[Chi]]) u\[Phi][\[Chi]] Y)/(l (1+l) \[Lambda] rp[\[Chi]] ut[\[Chi]])]],Function[\[Chi],Evaluate[(8 E^(-I m \[Phi]p[\[Chi]]) \[Pi] (-2 M+rp[\[Chi]]) (-((rp[\[Chi]]^2 ur[\[Chi]]^2)/ut[\[Chi]])+(-2 M+rp[\[Chi]])^2 ut[\[Chi]]) Y)/((1+\[Lambda]) rp[\[Chi]]^2 (3 M+\[Lambda] rp[\[Chi]]))]]];
		G=If[OddQ[l+m],Function[\[Chi],Evaluate[1/(l (1+l) \[Lambda] rp[\[Chi]] ut[\[Chi]]^2) 16  E^(-I m \[Phi]p[\[Chi]]) Sqrt[l+l^2-m (1+m)] \[Pi] u\[Phi][\[Chi]] (2 M ut[\[Chi]]^2+rp[\[Chi]] (ur[\[Chi]]^2-ut[\[Chi]]^2)-I m rp[\[Chi]]^2 ur[\[Chi]] u\[Phi][\[Chi]]) Y]],Function[\[Chi],Evaluate[-(1/((1+\[Lambda]) (3 M+\[Lambda] rp[\[Chi]])^2 ut[\[Chi]]))8 E^(-I m \[Phi]p[\[Chi]]) \[Pi] (-(((3 M^2+(1+\[Lambda]) rp[\[Chi]] (6 M+\[Lambda] rp[\[Chi]])) ur[\[Chi]]^2)/rp[\[Chi]])-2 I m (2 M-rp[\[Chi]]) (3 M+\[Lambda] rp[\[Chi]]) ur[\[Chi]] u\[Phi][\[Chi]]+1/(\[Lambda] rp[\[Chi]]^3) (-2 M+rp[\[Chi]]) (\[Lambda] (-2 M+rp[\[Chi]]) (15 M^2+\[Lambda] rp[\[Chi]] (6 M+(1+\[Lambda]) rp[\[Chi]])) ut[\[Chi]]^2+rp[\[Chi]]^3 (3 M+\[Lambda] rp[\[Chi]]) (M (3-3 m^2+5 \[Lambda])+\[Lambda] (-m^2+\[Lambda]) rp[\[Chi]]) u\[Phi][\[Chi]]^2)) Y]]];
		
		ZInToint = Function[\[Chi],Evaluate[(1-(2M)/rp[\[Chi]])^-1 PsiUp[\[Chi]]*G[\[Chi]]+((2M)/(rp[\[Chi]]-2M)^2 PsiUp[\[Chi]]-(1-(2M)/rp[\[Chi]])^-1 dPsiUp[\[Chi]])F[\[Chi]]]];
		ZUpToint = Function[\[Chi],Evaluate[(1-(2M)/rp[\[Chi]])^-1 PsiIn[\[Chi]]*G[\[Chi]]+((2M)/(rp[\[Chi]]-2M)^2 PsiIn[\[Chi]]-(1-(2M)/rp[\[Chi]])^-1 dPsiIn[\[Chi]])F[\[Chi]]]];
		
		(*t wrt \[Chi]*)
		dtd\[Chi] = Function[\[Chi], Evaluate[(M p^2)/((p-2-2e Cos[\[Chi]])(1+e Cos[\[Chi]])^2) (((p-2)^2-4e^2)/(p-6-2e Cos[\[Chi]]))^(1/2)]];
		
		(*Split into odd and even and half the integration period later once working?*)
		ZIn = \[CapitalOmega]r/(W 2\[Pi]) NIntegrate[ZInToint[\[Chi]]*Exp[I*\[Omega]*tp[\[Chi]]]dtd\[Chi][\[Chi]],{\[Chi],0,2\[Pi]}];
		ZUp = \[CapitalOmega]r/(W 2\[Pi]) NIntegrate[ZUpToint[\[Chi]]*Exp[I*\[Omega]*tp[\[Chi]]]dtd\[Chi][\[Chi]],{\[Chi],0,2\[Pi]}];

		<|"\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn|>
];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];

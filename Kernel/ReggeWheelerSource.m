(* ::Package:: *)

(* ::Title:: *)
(*ReggeWheelerSource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ReggeWheelerSource`", {"KerrGeodesics`KerrGeoOrbit`","KerrGeodesics`OrbitalFrequencies`","KerrGeodesics`FourVelocity`"}];


(* ::Subsection:: *)
(*Usage messages*)


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*ReggeWheelerPointParticleSource*)


(*Change from print to Error message*)
ReggeWheelerPointParticleSource[s_,l_,m_,n_, orbit_] :=
 Which[Abs[orbit["e"]] == 0 && orbit["Inclination"] == 1, Return[ReggeWheelerPointParticleSourceCircularCPMZM[s,l,m,orbit]], Abs[orbit["e"]] != 0&& orbit["Inclination"] == 1, Return[ReggeWheelerPointParticleSourceEccentricCPMZM[s,l,m,n,orbit]], True, Print["No point-particle source yet available for those parameters"]];


Format[ReggeWheelerSourceObject[assoc_]] := "ReggeWheelerSourceObject[<<>>]";

ReggeWheelerSourceObject[assoc_][string_] := assoc[string];


(* ::Subsection:: *)
(*Eccentric Orbits*)


(*the jump coditions are in the convolution module; currently developing for |s|=2 eccentric orbits using the CPM and ZM masterfunctions*)

ReggeWheelerPointParticleSourceEccentricCPMZM[s_,l_,m_,n_,orbit_]:=
	Module[{M=1,p,e,x,\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[Omega],assoc,traj},
		p=orbit["p"];
		e=orbit["e"];
		x=orbit["Inclination"];
		\[CapitalOmega]r=KerrGeoFrequencies[0,p,e,x]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)"];
		\[CapitalOmega]\[Phi]=KerrGeoFrequencies[0,p,e,x]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"];
		\[Omega]=m*\[CapitalOmega]\[Phi]+n*\[CapitalOmega]r;
		
		assoc=Association[
			"l"->l,
			"m"->m,
			"n"->n,
			"p"->p,
			"e"->e,
			"x"->x,
			"\[Omega]"->\[Omega],
			"\[CapitalOmega]r"->\[CapitalOmega]r,
			"\[CapitalOmega]\[Phi]"->\[CapitalOmega]\[Phi],
			"type"->"PointParticleEccentric"
		];
		ReggeWheelerSourceObject[assoc]
	];



(* ::Subsection::Closed:: *)
(*Circular Orbits*)


(*the jump coditions used in the convolution to find the particular solution
currently valid for s=2 circular orbits using the CPM and ZM masterfunctions*)

ReggeWheelerPointParticleSourceCircularCPMZM[s_,l_,m_,orbit_]:=
	If[EvenQ[l+m],
		ReggeWheelerPointParticleCircularZMEvenJump[l,m,orbit],
		ReggeWheelerPointParticleCircularCPMOddJump[l,m,orbit]
	];

ReggeWheelerPointParticleCircularCPMOddJump[l_,m_,orbit_]:=
	Module[{M=1,p,r0,\[CapitalOmega],ut,Y,F,G,\[Lambda]=((l+2)(l-1))/2,assoc},
		p=orbit["p"];
		r0=p;
		\[CapitalOmega]=KerrGeoFrequencies[0,p,0,1]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"];
		r0=p;
		Y=SphericalHarmonicY[l,m+1,\[Pi]/2,0];
		ut=Sqrt[r0/(r0-3M)];
		F=(16 Sqrt[(l-m) (1+l+m)] \[Pi] (-2 M+r0)^2 ut \[CapitalOmega] Y)/(l (1+l) r0 \[Lambda]);
		G=-((16 Sqrt[(l-m) (1+l+m)] \[Pi] (-2 M+r0) ut \[CapitalOmega] Y)/(l (1+l) r0 \[Lambda]));
		
		assoc=Association[
			"l"->l,
			"m"->m,
			"r0"->r0,
			"deltaPsi"->(F*r0^2)/(r0-2M)^2,
			"deltadPsidr"->(G*r0^2)/(r0-2M)^2+(F*r0*2M)/(r0-2M)^3,
			"type"->"PointParticleCircular"
		];
		ReggeWheelerSourceObject[assoc]
	];
	
ReggeWheelerPointParticleCircularZMEvenJump[l_,m_,orbit_]:=
	Module[{M=1,p,r0,\[CapitalOmega],ut,Y,F,G,\[Lambda]=((l+2)(l-1))/2,assoc},
		p=orbit["p"];
		r0=p;
		\[CapitalOmega]=KerrGeoFrequencies[0,p,0,1]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"];
		r0=p;
		Y=SphericalHarmonicY[l,m,\[Pi]/2,0];
		ut=Sqrt[r0/(r0-3M)];
		F=(8 \[Pi] (-2 M+r0)^3 ut Y)/(r0^2 (1+\[Lambda]) (3 M+r0 \[Lambda]));
		G=-(1/(r0^3 \[Lambda] (1+\[Lambda]) (3 M+r0 \[Lambda])^2))8\[Pi] (2 M-r0) ut (-r0^3 \[Lambda]^2 (1+\[Lambda])+M r0^2 \[Lambda]^2 (-4+m^2+\[Lambda])+2 M^2 r0 \[Lambda] (-9+3 m^2+2 \[Lambda])+3 M^3 (-3+3 m^2+5 \[Lambda])) Y;
		
		assoc=Association[
			"l"->l,
			"m"->m,
			"r0"->r0,
			"deltaPsi"->(F*r0^2)/(r0-2M)^2,
			"deltadPsidr"->(G*r0^2)/(r0-2M)^2+(F*r0*2M)/(r0-2M)^3,
			"type"->"PointParticleCircular"
		];
		ReggeWheelerSourceObject[assoc]
	];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];

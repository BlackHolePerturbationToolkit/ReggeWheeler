(* ::Package:: *)

(* ::Title:: *)
(*ReggeWheelerSource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`ReggeWheelerSource`"];


(* ::Subsection::Closed:: *)
(*Usage messages*)


ReggeWheelerSourceObject::usage = "ReggeWheelerSourceObject[assoc] an object which contains a Regge Wheeler source."
ReggeWheelerPointParticleSource::usage = "ReggeWheelerPointParticleSource[s, orbit] Point particle source for the Regge Wheeler equation."


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ReggeWheelerPointParticleSource*)


ReggeWheelerPointParticleSource[s_,l_,m_, orbit_] :=
 If[orbit["e"] == 0 && Abs[orbit["Inclination"]] == 1,
         Return[ReggeWheelerPointParticleSourceCircular[s,l,m,orbit]],
         Print["No point-particle source yet available for those parameters"];]


Format[ReggeWheelerSourceObject[assoc_]] := "ReggeWheelerSourceObject[<<>>]";

ReggeWheelerSourceObject[assoc_][string_] := assoc[string];


ReggeWheelerPointParticleCircularAKStressEnergy[l_,m_,r0_]:=
	Module[{E0,L0},
		L0=Sqrt[r0^2/(r0-3)];
		E0=Sqrt[(r0-2)(r0^2+L0^2)/r0^3];
		Association[
		"EA"->-16*Pi*(r0-2)*E0/r0^3*Conjugate[SphericalHarmonicY[l,m,Pi/2,0]],
		"EB"->0,
		"EC"->If[l>=1,-16*Pi*(r0-2)*L0/r0^4/(l(l+1))Conjugate[Derivative[0,0,1,0][SphericalHarmonicY][l,m,Pi/2,0]],0],
		"ED"->0,
		"EE"->0,
		"EF"->If[l>=2,-16*Pi*(r0-2)*L0^2/E0/r0^5/((l-1)*l*(l+1)*(l+2))*(2*Conjugate[Derivative[0,0,0,2][SphericalHarmonicY][l,m,Pi/2,0]]+l*(l+1)*Conjugate[SphericalHarmonicY[l,m,Pi/2,0]]),0],
		"EG"->0,
		"EH"->0,
		"EJ"->0,
		"EK"->0	
		]
	];


(*the jump coditions used in the convolution to find the particular solution
currently valid for s=2 circular orbits, but expandable for generic A-K s=2 point-particle stress-energies*)

ReggeWheelerPointParticleSourceCircular[s_,l_,m_,orbit_]:=
	If[EvenQ[l+m],
		ReggeWheelerPointParticleCircularAKEvenJump[l,m,orbit],
		ReggeWheelerPointParticleCircularAKOddJump[l,m,orbit]
	];

ReggeWheelerPointParticleCircularAKOddJump[l_,m_,orbit_]:=
	Module[{stressenergy,EC,EJ,r0,\[Omega],assoc},
		r0=orbit["p"];
		\[Omega]=m*Sqrt[1/r0^3];
		stressenergy=ReggeWheelerPointParticleCircularAKStressEnergy[l,m,r0];
		EC=stressenergy["EC"];
		EJ=stressenergy["EJ"];
		assoc=Association[
			"l"->l,
			"m"->m,
			"r0"->r0,
			"deltaPsi"->r0^3/(r0-2)*EC,
			"deltadPsidr"->-2*r0^2/(r0-2)^2*EC+r0^2/(r0-2)*EC-I*\[Omega]*r0^3/(r0-2)*EJ-3r0^2/(r0-2)*EC+r0^3/(r0-2)^2*EC,
			"type"->"PointParticleCircular"
		];
		ReggeWheelerSourceObject[assoc]
	];
	
ReggeWheelerPointParticleCircularAKEvenJump[l_,m_,orbit_]:=
	Module[{n,rm2M,np6M,term1,dterm1,coeff,term2,EA,ED,EF,EH,EK,stressenergy,r0,\[Omega],assoc},
		r0=orbit["p"];
		\[Omega]=m*Sqrt[1/r0^3];
		n=(l-1)(l+2);
		rm2M=r0-2;
		np6M=n*r0+6;
		stressenergy=ReggeWheelerPointParticleCircularAKStressEnergy[l,m,r0];
		EA=stressenergy["EA"];
		ED=stressenergy["ED"];
		EF=stressenergy["EF"];
		EH=stressenergy["EH"];
		EK=stressenergy["EK"];
		term1=-r0^4*EA/2/np6M/rm2M;
		dterm1=-r0^4*EA/2/np6M/rm2M*(4/r0-1/rm2M-n/np6M);
		coeff=2/r0/rm2M;
		term2=r0^3/4*(r0^2*n*(n-2)+r0*(14n-36)+96)*EA/(rm2M*np6M)^2+r0^4*I/2*\[Omega]*ED/np6M/rm2M+(n+2)r0^2/4*EF/rm2M-(n+2)r0^2/4*(2*EH+EK)/np6M;
		assoc=Association[
			"l"->l,
			"m"->m,
			"r0"->r0,
			"deltaPsi"->term1,
			"deltadPsidr"->-coeff*term1-dterm1+term2,
			"type"->"PointParticleCircular"
		];
		ReggeWheelerSourceObject[assoc]
	];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];

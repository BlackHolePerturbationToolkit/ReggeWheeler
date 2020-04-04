(* ::Package:: *)

BeginPackage["ReggeWheeler`ConvolveSource`"];

Begin["`Private`"];


ConvolveSource[RF_, SO_] :=Module[{s},
		s = RF["In"]["s"];
		
		If[Abs[s] == 2, 
			If[SO["type"]=="PointParticleCircular", Return[ConvolvePointParticleSourceCircular[s,RF,SO]], Print["Only circular orbit sources are implemented"]];
			,Print["Only the |s|=2 source is currently implemented"];
		];
	]


ConvolvePointParticleSourceCircular[(-2|2),RF_,SO_]:=
	Module[{l,m,r0,PsiIn,dPsiIn,PsiUp,dPsiUp,PsiOddIn,dPsiOddIndr,PsiOddUp,dPsiOddUpdr,Wronskian,deltadPsidr,deltaPsi,n,np6M,denom,conjdenom,b,c,\[Omega],rm2M,ZIn,ZUp,jump},
		l=SO["l"];
		m=SO["m"];
		r0=SO["r0"];
		\[Omega]=m*Sqrt[1/r0^3];
		rm2M=r0-2;
		If[EvenQ[l+m],
			PsiOddIn=RF["In"][r0];
			dPsiOddIndr=RF["In"]'[r0];
			PsiOddUp=RF["Up"][r0];
			dPsiOddUpdr=RF["Up"]'[r0];
			n=(l-1)*(l+2);
			np6M=n*r0+6;
			denom=n*(n+2)+12*I*\[Omega];
			conjdenom=n*(n+2)-12*I*\[Omega];
			b=((r0*n)^2*(np6M+2*r0)+12*r0^2*n+72*rm2M)/(r0^2*np6M);
			c=((r0*n)^2*(np6M+2*r0)+12*r0^2*n+36*rm2M)*12/r0/(r0*np6M)^2;
			PsiIn=(12*(1-2/r0)*dPsiOddIndr+b*PsiOddIn)/conjdenom;
			PsiUp=(12*(1-2/r0)*dPsiOddUpdr+b*PsiOddUp)/denom;
			dPsiIn=(-12*\[Omega]^2/(1-2/r0)*PsiOddIn+c*PsiOddIn+b*dPsiOddIndr)/conjdenom;
			dPsiUp=(-12*\[Omega]^2/(1-2/r0)*PsiOddUp+c*PsiOddUp+b*dPsiOddUpdr)/denom;
		,
			PsiIn=RF["In"][r0];
			dPsiIn=RF["In"]'[r0];
			PsiUp=RF["Up"][r0];
			dPsiUp=RF["Up"]'[r0];
		];
		
		Wronskian = PsiIn*dPsiUp - PsiUp*dPsiIn;
		deltaPsi=SO["deltaPsi"];
		deltadPsidr=SO["deltadPsidr"];
		ZIn = (PsiUp*deltadPsidr - deltaPsi*dPsiUp)/Wronskian;
		ZUp = (PsiIn*deltadPsidr - deltaPsi*dPsiIn)/Wronskian;
		<|"ZInf"->ZUp, "ZHor"->ZIn|>
	];


End[];
EndPackage[];

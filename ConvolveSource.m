(* ::Package:: *)

BeginPackage["ReggeWheeler`ConvolveSource`"];

Begin["`Private`"];


ConvolveSource[RF_, SO_] :=
	If[SO["type"]=="PointParticleCircular",ConvolvePointParticleSourceCircular[RF,SO],Nothing];


ConvolvePointParticleSourceCircular[RF_,SO_]:=
	Module[{l,m,r0,PsiIn,dPsiIn,PsiUp,dPsiUp,PsiOddIn,dPsiOddIndr,PsiOddUp,dPsiOddUpdr,Wronskian,deltadPsidr,deltaPsi,n,np6M,denom,conjdenom,b,c,\[Omega],rm2M,ZIn,ZUp,jump},
		l=SO["l"];
		m=SO["m"];
		r0=SO["r0"];
		\[Omega]=m*Sqrt[1/r0^3];
		rm2M=r0-2;
		If[EvenQ[l+m],
			PsiOddIn=RF["In"]["Psi"][r0];
			dPsiOddIndr=RF["In"]["dPsidr"][r0];
			PsiOddUp=RF["Up"]["Psi"][r0];
			dPsiOddUpdr=RF["Up"]["dPsidr"][r0];
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
			PsiIn=RF["In"]["Psi"][r0];
			dPsiIn=RF["In"]["dPsidr"][r0];
			PsiUp=RF["Up"]["Psi"][r0];
			dPsiUp=RF["Up"]["dPsidr"][r0];
		];
		Wronskian = PsiIn*dPsiUp - PsiUp*dPsiIn;
		deltaPsi=SO["deltaPsi"];
		deltadPsidr=SO["deltadPsidr"];
		ZIn = (PsiUp*deltadPsidr - deltaPsi*dPsiUp)/Wronskian;
		ZUp = (PsiIn*deltadPsidr - deltaPsi*dPsiIn)/Wronskian;
		<|"ZInf"->ZUp,"ZHor"->ZIn|>
	];


End[];
EndPackage[];

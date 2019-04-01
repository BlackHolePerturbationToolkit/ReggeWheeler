(* ::Package:: *)

(*use 'x' as r/M and \[Omega] to denote 'M\[Omega]'*)


BeginPackage["ReggeWheeler`NumericalIntegration`"];


Begin["`Private`"];


SetAttributes[PsiIn, {NumericFunction}];

PsiIn[s_, l_, \[Omega]_, xmin_, xmax_]:=
Module[{newassoc},
	If[\[Omega]==0,
		newassoc=<|
				"Psi"->Function[y,PsiInStaticOdd[s,l,y]],
				"dPsidr"->Function[y,Derivative[0,0,1][PsiInStaticOdd][s,l,y]],
				"xmin"->xmin,
				"xmax"->xmax
				|>;
		,
		If[s==2,
			Module[{boundaryconditions,soln,newAssoc,psiBC,dpsidxBC,xBC},
				boundaryconditions = ReggeWheelerInBC[s,l,\[Omega],$MachinePrecision];
				psiBC=boundaryconditions["Psi"];
				dpsidxBC=boundaryconditions["dPsidx"];
				xBC=boundaryconditions["xBC"];
				soln = Integrator[s,l,\[Omega],psiBC,dpsidxBC,xBC,xmax,ReggeWheelerPotential,$MachinePrecision];
				newassoc=<|
					"Psi"->Function[y,y1[y]/.soln],
					"dPsidr"->Function[y,y2[y]/.soln],
					"xmin"->xBC,
					"xmax"->xmax					
					|>;
			];,
			Nothing (*wait for further functionality*)
		];
	];
newassoc
];

(*Derivative[1][PsiIn[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][x_?InexactNumberQ] :=
 Module[{},
];

Derivative[n_Integer?Positive][PsiIn[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][r0_?InexactNumberQ] :=
 Module[{d2R, pderivs, R, r, i},
  d2R = -2/(x^2*(1-2/x))*Derivative[1][R][r] -(\[Omega]^2-ReggeWheelerPotential[s,l,x])/(1-2/x)^2*R[r];
  
  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> PsiIn[s, l, \[Omega]]'[r0],
    R[r] -> PsiIn[s, l, \[Omega]][r0], r -> r0}
];*)


SetAttributes[PsiUp, {NumericFunction}];

PsiUp[s_, l_, \[Omega]_, xmin_, xmax_]:=
	Module[{newassoc},
		If[\[Omega]==0,
			newassoc=<|
				"Psi"->Function[y,PsiUpStaticOdd[s,l,y]],
				"dPsidr"->Function[y,Derivative[0,0,1][PsiUpStaticOdd][s,l,y]],
				"xmin"->xmin,
				"xmax"->xmax
				|>;
			,
			If[s==2,
				Module[{boundaryconditions,soln,newAssoc,psiBC,dpsidxBC,xBC},
					boundaryconditions = ReggeWheelerUpBC[s,l,\[Omega],xmax,$MachinePrecision];
					psiBC=boundaryconditions["Psi"];
					dpsidxBC=boundaryconditions["dPsidx"];
					xBC=boundaryconditions["xBC"];
					soln = Integrator[s,l,\[Omega],psiBC,dpsidxBC,xBC,xmin,ReggeWheelerPotential,$MachinePrecision];
					newassoc= <|
						"Psi"->Function[y,y1[y]/.soln],
						"dPsidr"->Function[y,y2[y]/.soln],
						"xmin"->xmin,
						"xmax"->xBC
					|>
				],
			newassoc=Nothing (*wait for further functionality*)
			];
		];
		newassoc
	];

(*Derivative[1][PsiUp[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][x_?InexactNumberQ] :=
 Module[{},
];

Derivative[n_Integer?Positive][PsiUp[s_Integer, l_Integer, \[Omega]_, rmin_, rmax_]][r0_?InexactNumberQ] :=
 Module[{d2R, pderivs, R, r, i},
  d2R = -2/(x^2*(1-2/x))*Derivative[1][R][r] -(\[Omega]^2-ReggeWheelerPotential[s,l,x])/(1-2/x)^2*R[r];

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> PsiUp[s, l, \[Omega]]'[r0],
    R[r] -> PsiUp[s, l, \[Omega]][r0], r -> r0}
];*)


(*should this be in a module for y1 and y2?*)
Integrator[s_,l_,\[Omega]_,y1BC_,y2BC_,xBC_,xend_,potential_,precision_]:=
	NDSolve[
		{y1'[x]==y2[x],(1-2/x)^2*y2'[x]+2(1-2/x)/x^2*y2[x]+(\[Omega]^2-potential[s,l,x])*y1[x]==0,y1[xBC]==y1BC,y2[xBC]==y2BC},
		{y1,y2},
		If[xBC<xend,{x,xBC,xend},{x,xend,xBC}],
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		WorkingPrecision->precision,
		InterpolationOrder->All
	];


(*boundary conditions for odd-parity Regge-Wheeler eqn.*)
(*currently only working for s=2*)
(*SetAttributes[ReggeWheelerInBC, {NumericFunction}];
SetAttributes[ReggeWheelerUpBC, {NumericFunction}];*)

ReggeWheelerInBC[s_Integer,l_Integer,\[Omega]_,workingprecision_Integer]:=
	Module[{rm2M,p,ptrys,expeh,Dexpeh,done=False,delReh,Reh,rstar,drstardr,
	nmax,Xn,n,bk,denominator,f1=0,f2=0,f3=0,last,Bkm1=1,Bkm2=0,Bkm3=0,next,psi,dpsidr,om,precision=workingprecision+10,count},
		rm2M=4*^-2*\[Omega]^2/(l*(l+1));
		om=-\[Omega]; (* <-- this makes the sign of the exponential positive!*)
		p=0;
		ptrys=6;
		nmax=10;
		f1=0;
		f2=0;
		f3=0;
		Bkm1=1;
		Bkm2=0;
		Bkm3=0;
		expeh=1;
		Dexpeh=0;
		Xn=1;
		last=1*^33;
		count=0;
		While[!done && p<ptrys,
			delReh=rm2M;
			Reh=2+rm2M;
			p++;
			bk=1;
			For[n=1, n<=nmax, n++,
				Xn*=delReh;
				denominator=4(n+4*I*om)*n;
				f1=2(-12*I*om*(n-1)+l(l+1)-2n^2+5n-6)/denominator;
				f2=(-12*I*om*(n-2)+l(l+1)-n^2+5n-6)/denominator;
				f3=(-2(n-3)*I*om)/denominator;
				bk=N[f1*Bkm1+f2*Bkm2+f3*Bkm3,precision];
				Bkm3=Bkm2;
				Bkm2=Bkm1;
				Bkm1=bk;
				next=bk*Xn;
				If[(Abs[next]>Abs[last])&&(n>5),{count+=1,If[count>3,Break[]]}];
				If[n>6 && (RealExponent[N[(expeh+next)-expeh,precision]])<-(precision-5),{done=True,Break[]}];
				last=next;
				expeh+=next;
				Dexpeh+=n*next/delReh;
			];
			If[done==False,
				f1=0;
				f2=0;
				f3=0;
				Bkm1=1;
				Bkm2=0;
				Bkm3=0;
				expeh=1;
				Dexpeh=0;
				Xn=1;
				rm2M/=2;
				last=1*^33
				count=0;
			];
		];
	rstar = Reh+2*Log[rm2M/(2)];
	drstardr = Reh/rm2M;
	psi=Exp[I*om*rstar]*expeh;
	dpsidr=Exp[I*om*rstar]*Dexpeh+I*om*drstardr*psi;
	<|"Psi"->N[psi,precision],"dPsidx"->N[dpsidr,precision],"xBC"->N[Reh,precision]|>
]

ReggeWheelerUpBC[s_Integer,l_Integer,\[Omega]_,xmax_,workingprecision_Integer]:=
	Module[{An=1,Anm1=0,Anm2=0,Anm3=0,Nmax=75,NNmax=1000,n,nn,rstart,rstar,drstardr,increment=0,
	lastincrement=1*^40,S=0,lastS=0,dS=0,lastdS=0,count=0,r,rn,np,continue=True,precision=workingprecision+10,om,BCinc},
		(*rstarstart=xmax+2*Log[xmax/2-1]+5*Pi/Abs[om];*)
		rstart=xmax;
		om=-\[Omega];
		nn=1;
		r=rstart;
		rn=1;
		BCinc=10;
		While[continue&&nn<NNmax,
			For[n=0,n<=Nmax,n++,
				increment=An/rn;
				If[(Abs[increment]>Abs[lastincrement])&&(n>5),{count+=1,If[count>3,Break[]]}];
				lastS=S;
				lastdS=dS;
				S=S+increment;
				dS=dS-n*increment/r;
				If[n>4 && (RealExponent[N[S-lastS,precision]]<-(precision-5))&&(RealExponent[N[dS-(dS-n*increment/r),precision]]<-(precision-5)),{continue=False,Break[]}];
				np=n+1;
				rn=rn*r;
				lastincrement=Abs[increment];
				Anm3=Anm2;
				Anm2=Anm1;
				Anm1=An;
				An=N[-(4*np*(np-4)*Anm3+2*(l*(l+1)+3-(np-2)*(2np-1))*Anm2+((np-1)*np-l(l+1)-4*I*om*(np-1))*Anm1)/(2*I*om*np),precision];
			];
			If[continue==True,
			rstart+=BCinc;
			nn++;
			r=rstart;
			rn=1;
			An=1;
			Anm1=0;
			Anm2=0;
			Anm3=0;
			increment=0;
			lastincrement=1*^40;
			S=0;
			lastS=0;
			dS=0;
			lastdS=0;
			count=0;
			];
		];
	If[nn>=NNmax,{Print["The UP boundary condition loop ran out of steps!"],Abort[]}];
	rstar=rstart+2*Log[rstart/2-1];
	drstardr = rstart/(rstart-2);
	<|"Psi"->N[S*Exp[-I*om*rstart],precision],
	"dPsidx"->N[(-I*om*S*drstardr+dS)*Exp[-I*om*rstart],precision],
	"xBC"->N[rstart,precision]|>
]



(*radial potentials for ODE*)
SetAttributes[ReggeWheelerPotential,{NumericFunction}];
SetAttributes[ZerilliPotential,{NumericFunction}];

ReggeWheelerPotential[s_Integer,l_Integer,x_Numeric]:=(1-2/x)*(2(1-s^2)+l*(l+1)*x)/x^3;

(*only for spin s=2*)
ZerilliPotential[l_Integer,x_Numeric]:=
	Module[{n},
		n=(l-1)*(l+2)/2;
		2*(1-2/x)*(9+9*n*x+n^2*x^2*(3+x)+n^3*x^3)/(x^3(3+n*x)^2)
	];


(*analytic solutions to the homogeneous static ReggeWheeler Eq. (see e.g. Field, Hesthaven, and Lau, PRD81, 124030 (2010)*)
SetAttributes[PsiUpStaticOdd,{NumericFunction}]
SetAttributes[PsiInStaticOdd,{NumericFunction}]

PsiUpStaticOdd[s_Integer,l_Integer,x_]:=(x/2)^(-l)*Hypergeometric2F1[l+s+1,l-s+1,2*(l+1),2/x];

PsiInStaticOdd[s_Integer,l_Integer,x_]:=(x/2)^(-l)*Hypergeometric2F1[l+s+1,l-s+1,1,(x-2)/x];

(*SetAttributes[PsiUpStaticEven,{NumericFunction}]
SetAttributes[PsiInStaticEven,{NumericFunction}]*)

(*even-parity static solutions for s=2 given by the intertwining operatior*)
PsiUpStaticEven[s_Integer,l_Integer,x_]:=
	If[s==2,(*handle spin-2 case*)
		Module[{n},
			n=(l-1)*(l+2)/2;
			2*(1-2/x)*Derivative[0,0,1][PsiUpStaticOdd][s,l,x]+(2*n*(n+1)/3+6*(x-2)/(x^2*(3+n*x)))*PsiUpStaticOdd[s,l,x]
		]
	,
		Nothing
	];
	
PsiInStaticEven[s_Integer,l_Integer,x_]:=
	If[s==2,(*handle spin-2 case*)
		Module[{n},
			n=(l-1)*(l+2)/2;
			2*(1-2/x)*Derivative[0,0,1][PsiInStaticOdd][s,l,x]+(2*n*(n+1)/3+6*(x-2)/(x^2*(3+n*x)))*PsiInStaticOdd[s,l,x]
		]
	,
		Nothing
	];


(*useful functions?*)
rfromrstar=Function[rstar,2 (1+ProductLog[Sqrt[E^(-2+rstar)]])];


End[]
EndPackage[];

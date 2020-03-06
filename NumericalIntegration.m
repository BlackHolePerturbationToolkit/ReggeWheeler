(* ::Package:: *)

(*use 'x' as r/M and \[Omega] to denote 'M\[Omega]'*)


BeginPackage["ReggeWheeler`NumericalIntegration`"];


Begin["`Private`"];


SetAttributes[Psi, {NumericFunction}];

Psi[s_, l_, (0|0.), bc_][{xmin_, xmax_}] :=
 If[bc == "In",
   Function[y, PsiInStaticOdd[s,l,y]],
   Function[y, PsiUpStaticOdd[s,l,y]]
 ];

Psi[s_, l_, \[Omega]_, "In"][xmax_?NumericQ] := Psi[s, l, \[Omega], "In"][{Automatic, xmax}];
Psi[s_, l_, \[Omega]_, "Up"][xmin_?NumericQ] := Psi[s, l, \[Omega], "Up"][{xmin, Automatic}];

Psi[s_, l_, \[Omega]_, bc_][{xmin_, xmax_}] :=
 Module[{bcFunc, psiBC, dpsidxBC, xBC, xMin, xMax, soln},
  If[s==2,
    bcFunc = Lookup[<|"In" -> ReggeWheelerInBC, "Up" -> ReggeWheelerUpBC|>, bc];
    {psiBC, dpsidxBC, xBC} = bcFunc[s, l, \[Omega], $MachinePrecision];
    If[bc === "In" && xmin === Automatic, xMin = xBC, xMin = xmin];
    If[bc === "Up" && xmax === Automatic, xMax = xBC, xMax = xmax];
    soln = Integrator[s, l, \[Omega], psiBC, dpsidxBC, xBC, xMin, xMax, ReggeWheelerPotential, $MachinePrecision];
  ,
    soln = Function[{x}, $Failed] (*wait for further functionality*)
  ];
  soln
];

Psi[s_, l_, \[Omega]_, bc_][All] :=
 Module[{bcFunc, psiBC, dpsidxBC, xBC, xMin, xMax, soln},
  If[s==2,
    bcFunc = Lookup[<|"In" -> ReggeWheelerInBC, "Up" -> ReggeWheelerUpBC|>, bc];
    {psiBC, dpsidxBC, xBC} = bcFunc[s, l, \[Omega], $MachinePrecision];
    soln = Function[{x}, Evaluate[Integrator[s, l, \[Omega], psiBC, dpsidxBC, xBC, Min[x, xBC], Max[x, xBC], ReggeWheelerPotential, $MachinePrecision][x]]];
  ,
    soln = Function[{x}, $Failed] (*wait for further functionality*)
  ];
  soln
]


(*should this be in a module for y1 and y2?*)
Integrator[s_,l_,\[Omega]_,y1BC_,y2BC_,xBC_,xmin_?NumericQ,xmax_?NumericQ,potential_,precision_]:=Module[{y1,y2,x},
	NDSolveValue[
		{y1'[x]==y2[x],(1-2/x)^2*y2'[x]+2(1-2/x)/x^2*y2[x]+(\[Omega]^2-potential[s,l,x])*y1[x]==0,y1[xBC]==y1BC,y2[xBC]==y2BC},
		y1,
		{x, xmin, xmax},
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		WorkingPrecision->precision,
		InterpolationOrder->All
		]
	]


(*boundary conditions for odd-parity Regge-Wheeler eqn.*)
(*currently only working for s=2*)
(*SetAttributes[ReggeWheelerInBC, {NumericFunction}];
SetAttributes[ReggeWheelerUpBC, {NumericFunction}];*)

ReggeWheelerInBC[s_Integer,l_Integer,\[Omega]_,workingprecision_]:=
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
				bk=-(((-3+n) \[Omega] Bkm3)/(2 n (I n+4\[Omega])))+(I (l+l^2-(-2+n) (-3+n-12 I \[Omega])) Bkm2)/(4 n (I n+4 \[Omega]))-((2-l-l^2+2 n^2+s^2+12 I \[Omega]+n (-5-12 I M \[Omega])) Bkm1)/(2 n (n-4 I \[Omega]));
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
				last=1*^33;
				count=0;
			];
		];
	rstar = Reh+2*Log[rm2M/(2)];
	drstardr = Reh/rm2M;
	psi=Exp[I*om*rstar]*expeh;
	dpsidr=Exp[I*om*rstar]*Dexpeh+I*om*drstardr*psi;
	{N[psi,precision], N[dpsidr,precision], N[Reh,precision]}
]

ReggeWheelerUpBC[s_Integer,l_Integer,\[Omega]_,workingprecision_]:=
	Module[{An=1,Anm1=0,Anm2=0,Nmax=75,NNmax=1000,n,nn,rstart,rstar,drstardr,increment=0,
	lastincrement=1*^40,S=0,lastS=0,dS=0,lastdS=0,count=0,r,rn,np,continue=True,precision=workingprecision+10,om,BCinc},
		(*rstarstart=xmax+2*Log[xmax/2-1]+5*Pi/Abs[om];*)
		rstart=10000;
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
				Anm2=Anm1;
				Anm1=An;
				An=N[(I(-1+np-s) (-1+np+s) Anm2)/(np \[Omega])+(I (l+l^2+np-np^2) Anm1)/(2 np \[Omega]),precision];
			];
			If[continue==True,
			rstart+=BCinc;
			nn++;
			r=rstart;
			rn=1;
			An=1;
			Anm1=0;
			Anm2=0;
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
	{N[S*Exp[-I*om*rstar],precision],
	 N[(-I*om*S*drstardr+dS)*Exp[-I*om*rstar],precision],
	 N[rstart,precision]}
]



(*radial potentials for ODE*)
(*SetAttributes[ReggeWheelerPotential,{NumericFunction}];*)
SetAttributes[ZerilliPotential,{NumericFunction}];

ReggeWheelerPotential[s_Integer,l_Integer,x_]:=(1-2/x)*(2(1-s^2)+l*(l+1)*x)/x^3;

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

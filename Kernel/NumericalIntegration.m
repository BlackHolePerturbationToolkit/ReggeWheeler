(* ::Package:: *)

(* ::Title:: *)
(*NumericalIntegration*)


(* ::Section::Closed:: *)
(*Create Package*)


(*use 'x' as r/M and \[Omega] to denote 'M\[Omega]'*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["ReggeWheeler`NumericalIntegration`"];


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section:: *)
(*Radial solutions*)


(* ::Subsection::Closed:: *)
(*Radiative modes*)


SetAttributes[Psi, {NumericFunction}];

Psi[s_, l_, \[Omega]_, "In", ndsolveopts___][xmax_?NumericQ] := Psi[s, l, \[Omega], "In", ndsolveopts][{Automatic, xmax}];
Psi[s_, l_, \[Omega]_, "Up", ndsolveopts___][xmin_?NumericQ] := Psi[s, l, \[Omega], "Up", ndsolveopts][{xmin, Automatic}];

Psi[s_, l_, \[Omega]_, bc_, ndsolveopts___][{xmin_, xmax_}] :=
 Module[{bcFunc, psiBC, dpsidxBC, xBC, xMin, xMax, soln},
    bcFunc = Lookup[<|"In" -> ReggeWheelerInBC, "Up" -> ReggeWheelerUpBC|>, bc];
    {psiBC, dpsidxBC, xBC} = bcFunc[s, l, \[Omega], Lookup[{ndsolveopts}, WorkingPrecision, Precision[\[Omega]]]];
    If[bc === "In" && xmin === Automatic, xMin = xBC, xMin = xmin];
    If[bc === "Up" && xmax === Automatic, xMax = xBC, xMax = xmax];
    soln = Integrator[s, l, \[Omega], psiBC, dpsidxBC, xBC, xMin, xMax, ReggeWheelerPotential, ndsolveopts]
];

Psi[s_, l_, \[Omega]_, bc_, ndsolveopts___][All] :=
 Module[{bcFunc, psiBC, dpsidxBC, xBC, xMin, xMax, soln},
    bcFunc = Lookup[<|"In" -> ReggeWheelerInBC, "Up" -> ReggeWheelerUpBC|>, bc];
    {psiBC, dpsidxBC, xBC} = bcFunc[s, l, \[Omega], Lookup[{ndsolveopts}, WorkingPrecision, Precision[\[Omega]]]];
    soln = AllIntegrator[s, l, \[Omega], psiBC, dpsidxBC, xBC, ReggeWheelerPotential, ndsolveopts]
];

Psi[s_, l_, \[Omega]_, bc_, ndsolveopts___][None] := $Failed;


(*should this be in a module for y1 and y2?*)
Integrator[s_,l_,\[Omega]_,y1BC_,y2BC_,xBC_,xmin_?NumericQ,xmax_?NumericQ,potential_,ndsolveopts___]:=Module[{y1,y2,x},
	Quiet[NDSolveValue[
		{y1'[x]==y2[x],(1-2/x)^2*y2'[x]+2(1-2/x)/x^2*y2[x]+(\[Omega]^2-potential[s,l,x])*y1[x]==0,y1[xBC]==y1BC,y2[xBC]==y2BC},
		y1,
		{x, xmin, xmax},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	]


(*should this be in a module for y1 and y2?*)
AllIntegrator[s_,l_,\[Omega]_,y1BC_,y2BC_,xBC_,potential_,ndsolveopts___][xval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,x},
	Quiet[NDSolveValue[
		{y1'[x]==y2[x],(1-2/x)^2*y2'[x]+2(1-2/x)/x^2*y2[x]+(\[Omega]^2-potential[s,l,x])*y1[x]==0,y1[xBC]==y1BC,y2[xBC]==y2BC},
		y1[xval],
		{x, Min[xBC,xval], Max[xBC,xval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];

Derivative[n_][AllIntegrator[s_,l_,\[Omega]_,y1BC_,y2BC_,xBC_,potential_,ndsolveopts___]][xval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,x},
	Quiet[NDSolveValue[
		{y1'[x]==y2[x],(1-2/x)^2*y2'[x]+2(1-2/x)/x^2*y2[x]+(\[Omega]^2-potential[s,l,x])*y1[x]==0,y1[xBC]==y1BC,y2[xBC]==y2BC},
		Derivative[n][y1][xval],
		{x, Min[xBC,xval], Max[xBC,xval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


(* ::Subsection:: *)
(*Boundary Conditions*)


(*boundary conditions for odd-parity Regge-Wheeler eqn.*)
(*currently only working for s=2*)
(*SetAttributes[ReggeWheelerInBC, {NumericFunction}];
SetAttributes[ReggeWheelerUpBC, {NumericFunction}];*)

ReggeWheelerMSTBC[s1_Integer, l1_Integer, \[Omega]1_, workingprecision_]:=
	Block[{s=s1,l=l1,\[Omega]=\[Omega]1,resIN,dresIN,resUP,dresUP,r,rin,rout, RW},
	
	rin = 2 + 10^-5;
	rout = 100\[Omega]^-1;
	
	RW = ReggeWheelerRadial[s,l,\[Omega],PrecisionGoal->workingprecision, AccuracyGoal->Infinity,Method->{"MST"}];
	
	resIN = RW["In"][rin];
	dresIN = D[RW["In"][r],r]/.r->rin;
	resUP = RW["Up"][rout];
	dresUP = D[RW["Up"][r],r]/.r->rout;
	
	<|"UpBC"->{resUP,dresUP,rout}, "InBC"->{resIN,dresIN,rin}|>
	
]

f=1-2/r;
RW[U_]:=f^2 D[\[Psi][r],{r,2}]+f D[f,r] D[\[Psi][r],r]+(\[Omega]^2-U)\[Psi][r]
Uodd[r_]:=f/r^2 (l(l+1)+(2M(1-s^2))/r);

rs[r_]:=r+2 Log[r/2-1];

ReggeWheelerInBC[s1_Integer, l1_Integer, \[Omega]1_, workingprecision_]:=
	Block[{s=s1,l=l1,i,A,a,n,\[Omega]=\[Omega]1,res,M=1,fH,err,r,rin=2+10^-5},
A[n_]:=((-1+2 n-n^2+s^2) a[-2+n]+(1+l+l^2-2 n+2 n^2-s^2) a[-1+n])/(n (n-4 I M \[Omega]));

a[-2]=0;
a[-1]=0;
a[0]=1;

fH=1-(2M)/r;

err=1;

res = Exp[-I \[Omega] rs[r]];
i=1;
While[err > 10^-workingprecision,

a[i]=A[i];

res+=Exp[-I \[Omega] rs[r]]a[i]fH^i;

err=Abs[RW[Uodd[r]]]/.{\[Psi][r]->res,Derivative[a_][\[Psi]][r]:>D[res,{r,a}]}/.r->rin;

i++;

];

{res,D[res,r],rin}/.r->rin
]

ReggeWheelerUpBC[s1_Integer,l1_Integer,\[Omega]1_,workingprecision_]:= Block[{s=s1,l=l1,A,a,n,\[Omega]=\[Omega]1,res,M=1,err,r,dres,d2res,rout, rsout,i},

(*This ensures the boundary is placed in the wavezone*)
rout =100\[Omega]^-1;

A[n_]:=(I (2 M \[Omega](1-2 n+n^2-s^2) a[-2+n]+(l+l^2+n-n^2) a[-1+n]))/(2 n );

a[-2]=0;
a[-1]=0;
a[0]=1;

err=1;

rsout = rs[rout];

res = Exp[I \[Omega] rsout];
dres = -((I E^(I \[Omega] rsout) rout \[Omega])/(2 M-rout));
d2res=-((E^(I \[Omega] rsout) \[Omega] (2 I M+rout^2 \[Omega]))/(-2 M+rout)^2);
i=1;
While[err > 10^-workingprecision,

a[i]=A[i];
res+=Exp[I \[Omega] rsout] a[i]/( \[Omega] rout)^i;
dres+=(E^(I \[Omega] rsout) (rout \[Omega])^-i (2 i M-i rout+I rout^2 \[Omega]) a[i])/(rout (-2 M+rout));
d2res+= (E^(I \[Omega] rsout) (rout \[Omega])^-i (i^2 (-2 M+rout)^2-rout^2 \[Omega] (2 I M+rout^2 \[Omega])+i (2 M-rout) (2 M+rout (-1+2 I rout \[Omega]))) a[i])/(rout^2 (-2 M+rout)^2);

err=Abs[RW[Uodd[r]]]/.{\[Psi][r]->res,\[Psi]'[r]->dres,\[Psi]''[r]->d2res}/.r->rout;

i++;

(*The infinity BCs is an asymptotic series so it might not converge.*)
If[i > 100, Break[]];

];

{res,dres,rout}
]


(* ::Subsection::Closed:: *)
(*Radial Potentials*)


(*SetAttributes[ReggeWheelerPotential,{NumericFunction}];*)
SetAttributes[ZerilliPotential,{NumericFunction}];

ReggeWheelerPotential[s_Integer,l_Integer,x_]:=(1-2/x)(2(1-s^2)+l(l+1)x)/x^3;

(*only for spin s=2*)
ZerilliPotential[l_Integer,x_Numeric]:=
	Module[{n},
		n=(l-1)*(l+2)/2;
		2*(1-2/x)*(9+9*n*x+n^2*x^2*(3+x)+n^3*x^3)/(x^3(3+n*x)^2)
	];


(* ::Subsection::Closed:: *)
(*Useful functions?*)


rfromrstar=Function[rstar,2 (1+ProductLog[Sqrt[E^(-2+rstar)]])];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];

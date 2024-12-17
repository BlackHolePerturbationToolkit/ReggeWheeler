(* ::Package:: *)

(* ::Chapter:: *)
(*ReggeWheelerHyperboloidal*)


(* ::Subsubsection:: *)
(*For details on solving the Regge-Wheeler-Zerilli equations using hyperboloidal compactification, see arXiv : gr - qc/2202.01794 and arXiv : gr - qc/2411.14976 .*)


(* ::Section:: *)
(*Create Package*)


(* ::Subsection:: *)
(*Begin Package*)


BeginPackage["ReggeWheeler`Hyperboloidal`",
  {"ReggeWheeler`ReggeWheelerSource`",
   "ReggeWheeler`ReggeWheelerRadial`",
   "ReggeWheeler`ConvolveSource`",
   "KerrGeodesics`KerrGeoOrbit`",
   "KerrGeodesics`OrbitalFrequencies`",
   "SpinWeightedSpheroidalHarmonics`"}
];


(* ::Subsection:: *)
(*Begin Private Section*)


Begin["`Private`"];


(* ::Section:: *)
(*Coordinate Transformation*)


DomainMapping[r0_,x_,X_]:= 
		Module[{M, r, \[Sigma]p, \[Sigma]grid1, \[Sigma]grid2, AB1, AB2, map1, map2, InvMap1, InvMap2, a, b},
			
			(* Initial setup *)
			M =1;
			\[Sigma]p=(2M)/r0;
			
			\[Sigma]grid1 = {0,\[Sigma]p};
			\[Sigma]grid2 = {\[Sigma]p,1};
			
			(* Solving for the transform coefficients *);
			AB1 =Solve[{a*\[Sigma]grid1[[1]]+b==-1,a*\[Sigma]grid1[[2]]+b==1}][[1]];
			AB2 = Solve[{a*\[Sigma]grid2[[1]]+b==-1,a*\[Sigma]grid2[[2]]+b==1}][[1]];
			
			(* Map to hyp. coords. *);
			map1= a*x+b/.AB1;
			map2=a*x+b/.AB2;
			
			(* Map to Sch. coords *);
			InvMap1 = (X-b)/a/.AB1;
			InvMap2=(X-b)/a/.AB2;
			
		Return[{map1,map2,InvMap1,InvMap2}]
];


(* ::Section:: *)
(*Utility functions*)


(* Schwarzschild f(r) *)
	f[r_,M_]:= 1-2 M/r;
	
(* Fns. defined for ease of repeated use *)
	L[l_] := l(l+1); 
	\[Lambda][l_]:=((l+2)(l-1))/2; 
	\[CapitalLambda][r_,l_,M_]:=\[Lambda][l] + (3M)/r;
	p[r_]:=(8*\[Pi])/r^2;
	q[r_,l_,M_]:= f[r,M]^2/((\[Lambda][l]+1)\[CapitalLambda][r,l,M]);
	Y\[Phi][l_,m_,\[Theta]_]:= Conjugate[(((D[SphericalHarmonicY[l,m,\[Theta],\[Phi]],{\[Phi],2}])/.\[Phi]->0)+Sin[\[Theta]]Cos[\[Theta]]((D[SphericalHarmonicY[l,m,x,0],x])/.x->\[Theta])
					+L[l]/2 Sin[\[Theta]]^2 SphericalHarmonicY[l,m,\[Theta],0])];
	X\[Theta][l_,m_,\[Theta]_]:= Conjugate[Sin[\[Theta]]D[SphericalHarmonicY[l,m,x,0],x]/.x->\[Theta]];
	
(* azimuthal orbital frequency *)
	\[CapitalOmega][r_,M_] := 1/r Sqrt[M/r]; 
	
(* angular frequency *)
	\[Omega][\[CapitalOmega]_,m_]:= m \[CapitalOmega];
	
(* orbital angular momentum per unit mass *)
	L0[r_,M_]:=r Sqrt[M/(r-3M)]; 
	
(* orbital energy per unit mass *)
	\[ScriptCapitalE][r_,M_]:=(r-2M)/Sqrt[r(r-3M)];
	
(* Even/Odd parity potentials *)	
	VeffEven[r_,l_,M_]:=f[r,M]/(r^2 \[CapitalLambda][r,l,M]^2) (2\[Lambda][l]^2 (\[Lambda][l]+1+(3M)/r)+18 M^2/r^2 (\[Lambda][l]+M/r));
	VeffOdd[r_,l_,M_]:= f[r,M]/r^2 (L[l] - 6 M/r); 

(* RWZ gauge source terms *)
	SourceTerm1Even[r_,M_,l_,m_,\[Theta]_]:=(( p[r]q[r,l,M] \[ScriptCapitalE][r,M])/(r f[r,M]\[CapitalLambda][r,l,M]) (L0[r,M]^2/\[ScriptCapitalE][r,M]^2 f[r,M]^2 \[CapitalLambda][r,l,M]-(\[Lambda][l](\[Lambda][l]+1)r^2+6\[Lambda][l]M r +15M^2))
									Conjugate[SphericalHarmonicY[l,m,\[Theta],0]]-(4p[r]L0[r,M]^2 f[r,M]^2)/(r \[ScriptCapitalE][r,M]) (l-2)!/(l+2)! Y\[Phi][l,m,\[Theta]]);
	SourceTerm2Even[r_,M_,l_,m_,\[Theta]_]:=(p[r] q[r,l,M]r^2 \[ScriptCapitalE][r,M])Conjugate[SphericalHarmonicY[l,m,\[Theta],0]];
	SourceTerm1Odd[r_,M_,l_,m_,\[Theta]_]:= -((2 p[r]f[r,M]L0[r,M])/(\[Lambda][l]*L[l]))X\[Theta][l,m,\[Theta]];
	SourceTerm2Odd[r_,M_,l_,m_,\[Theta]_]:= (2 p[r] r f[r,M]^2 L0[r,M])/(\[Lambda][l]*L[l])X\[Theta][l,m,\[Theta]];

(* Hyperboloidal height fn. *)
	H[\[Sigma]_] :=1/2 (Log[1-\[Sigma]]-1/\[Sigma]+Log[\[Sigma]]); 
	
(* Hyperboloidal scaling factor *)
	Z[\[Sigma]_,\[Xi]_]:= 1/2 E^(\[Xi] H[\[Sigma]]);
	
(* Hyperboloidal re-scaling factor *)
	\[ScriptCapitalF][r_,M_,\[Xi]_]:= f[r,M] 1/(2r^2) E^(\[Xi] H[2/r]);
	
(* Coefficients of hyperboloidal master eqn. operator *)
	\[Alpha]2[\[Sigma]_] := \[Sigma]^2 (1-\[Sigma]);
	\[Alpha]1[\[Sigma]_,\[Xi]_]:= \[Sigma](2-3\[Sigma])+\[Xi](1-2\[Sigma]^2);
	(* Distinct coefficients for the even/odd parity equations *)
	\[Alpha]0Odd[\[Sigma]_,l_,M_,\[Xi]_]:= -(\[Xi]^2 (1+\[Sigma])+2\[Xi] \[Sigma]+(4M^2 VeffOdd[2/\[Sigma],l,M])/((1-\[Sigma])\[Sigma]^2)); 
	\[Alpha]0Even[\[Sigma]_,l_,M_,\[Xi]_]:= -(\[Xi]^2 (1+\[Sigma])+2\[Xi] \[Sigma]+(4M^2 VeffEven[2/\[Sigma],l,M])/((1-\[Sigma])\[Sigma]^2));


(* ::Section:: *)
(*Spectral ODE Solver*)


Options[HyperboloidalSolver]={"GridPoints"->32};


HyperboloidalSolver[r0_, l_, m_, Xgrid_, opts:OptionsPattern[]]:=Module[
	{npts, M, \[Theta], \[Xi], \[Sigma]p, prec, map1, map2, InvMap1, InvMap2, x, X, \[Phi], S1, S2, ansatz, Dansatz, D2ansatz, Dmap1, Dmap2, 
	cs1, cs2, cs, DH, A, B, ansatz1D1, ansatz1D2, ansatz2D1, ansatz2D2, \[Alpha]21, 
	\[Alpha]22, \[Alpha]11, \[Alpha]12, \[Alpha]01, \[Alpha]02, BCs, BCsRHS, D\[Alpha]2, ODEs, juncs, juncs1,
	juncs2, fill, juncsRHS, dom1, dom2, Mat, LARHS, 
	sols, sols2, csNew, sol1, sol2, map1New, map2New, y, sol1New, sol2New, poly1, poly2},
	(*Module internal number of points on grid (maybe unnecessary)*)
		npts = OptionValue["GridPoints"];
		
	(* Initial setup *)
		M = 1;
		\[Theta] = \[Pi]/2;
		\[Xi] = -I \[Omega][\[CapitalOmega][r0,M],m] 4M;
		(* Source radial position in hyp coords *)
		\[Sigma]p = 2/r0;
		(* Setting working precision *)
		prec = Precision[r0];
		
	(* Coordinate mappings based off source position *)
		{map1, map2, InvMap1, InvMap2} = DomainMapping[r0,x,X];
		
	(* Source terms *)
		{S1,S2} = If[EvenQ[l+m],
						{SourceTerm1Even[r0,M,l,m,\[Theta]],SourceTerm2Even[r0,M,l,m,\[Theta]]},
						{SourceTerm1Odd[r0,M,l,m,\[Theta]],SourceTerm2Odd[r0,M,l,m,\[Theta]]}];
		
	(* Chebyshev polynomial as ansatz, using number of points on grid to determine length *)
		ansatz = Table[ChebyshevT[i,X],{i,0,Length[Xgrid]+1}];
		
	(* Derivatives of ansatz *)
		{Dansatz,D2ansatz,Dmap1,Dmap2} = {D[ansatz,X],D[ansatz,{X,2}],D[map1,x],D[map2,x]};
		
	(* Initialising table of weight coefficients *)
		{cs1,cs2} = {Table[Subscript[c, i, 1],{i,0,Length[Xgrid]+1}],Table[Subscript[c, i, 2],{i,0,Length[Xgrid]+1}]};
		cs = Join[cs1,cs2];
	
	(*Derivative of height fn for convenience*)
		DH =(D[H[x],x])/.x->\[Sigma]p;
		
	(* Hyperboloidal junction condition coefficients *)
		{A,B} = { (2E^(-\[Xi] H[\[Sigma]p]))/(1-\[Sigma]p) (2 S1+\[Sigma]p^2/(1-\[Sigma]p) (1-\[Xi](1-\[Sigma]p)(DH))S2),-((2 \[Sigma]p^2 E^(-\[Xi] H[\[Sigma]p]))/(1-\[Sigma]p))S2};
		
	(* Transforming ansatz derivatives  *)	
		{ansatz1D1, ansatz2D1, ansatz1D2, ansatz2D2} = {Dmap1*Dansatz, Dmap2*Dansatz, Dmap1^2*D2ansatz, Dmap2^2*D2ansatz};
		
	(* Obtaining master fn operator coefficients on Chebyshev-Gauss-Lobatto grid *)	
		{\[Alpha]21, \[Alpha]22, \[Alpha]11, \[Alpha]12, \[Alpha]01, \[Alpha]02} = If[EvenQ[l+m],
										{\[Alpha]2[InvMap1], \[Alpha]2[InvMap2], \[Alpha]1[InvMap1,\[Xi]], \[Alpha]1[InvMap2,\[Xi]], \[Alpha]0Even[InvMap1,l,M,\[Xi]], \[Alpha]0Even[InvMap2,l,M,\[Xi]]},
										{\[Alpha]2[InvMap1], \[Alpha]2[InvMap2], \[Alpha]1[InvMap1,\[Xi]], \[Alpha]1[InvMap2,\[Xi]], \[Alpha]0Odd[InvMap1,l,M,\[Xi]], \[Alpha]0Odd[InvMap2,l,M,\[Xi]]}];
		
	(* Defining boundary conditions, junction conditions and the ODE at every other grid point *)
		{BCs, ODEs, juncs1} = {{((\[Alpha]11 ansatz1D1 + \[Alpha]01 ansatz)/.X -> -1),((\[Alpha]12 ansatz2D1 + \[Alpha]02 ansatz)/.X -> 1)},
								{((\[Alpha]21 ansatz1D2 + \[Alpha]11 ansatz1D1 + \[Alpha]01 ansatz)/.X -> Xgrid),((\[Alpha]22 ansatz2D2 + \[Alpha]12 ansatz2D1 + \[Alpha]02 ansatz)/.X -> Xgrid)},
								{((cs2 . ansatz/.X -> -1) - (cs1 . ansatz/.X -> 1)),((cs2 . ansatz2D1/.X -> -1) - (cs1 . ansatz1D1/.X -> 1))}};
								
	 (* Obtaining values for junction conditions *)
		juncs2 = Table[Coefficient[juncs1[[j]],cs[[i]]], {j, Length[juncs1]}, {i, Length[cs]}];
	
	(* Placeholder column of zeroes for matrix construction (replaced by junction condition vectors in full matrix ) *)
		fill = ConstantArray[0, Length[ansatz]];
		
	(* Constructing matrices to be inverted for each domain *);
		dom1 = MapThread[Append, {MapThread[Prepend,{ODEs[[1]],BCs[[1]]}],fill}];
		dom2 = MapThread[Append,{MapThread[Prepend,{ODEs[[2]],fill}],BCs[[2]]}];
		
	(* Constructing overall block diagonal matrix to be inverted *)
		Mat = If[$VersionNumber >14.0, BlockDiagonalMatrix[{dom1,dom2}], ArrayFlatten[{{dom1,0},{0,dom2}}]];
		(*Setting junction conditions*)
		{Mat[[;;,Length[Mat]/2]], Mat[[;;,Length[Mat]/2+1]]} = juncs2;

	(*Initialising values for RHS of OVERALL equation Mx = b.*);
		BCsRHS = {0,0};
		D\[Alpha]2 = (D[\[Alpha]2[x],x])/.x->\[Sigma]p;
		juncsRHS = {B/\[Alpha]2[\[Sigma]p ],(A/\[Alpha]2[\[Sigma]p ]+B (D\[Alpha]2 -\[Alpha]1[\[Sigma]p,\[Xi] ])/\[Alpha]2[\[Sigma]p ]^2)};
		
	(* Constructing RHS vector *);
		LARHS = Join[{BCsRHS[[1]]},ConstantArray[0,Length[Xgrid]],juncsRHS,ConstantArray[0,Length[Xgrid]],{BCsRHS[[2]]}];
		
	(*Inverting matrix to obtain weight coefficients *);
		sols = LARHS . SetPrecision[Inverse[Mat],prec];
	
	
	(* Creating table of replacement rules for weight coefficients *);
		csNew = Table[cs[[i]]->sols[[i]],{i,1,Length[cs]}];
		
	(* Creating final polynomial *);
		sol1[x_]:= (cs1 . ansatz/.csNew)/.X->x;
		sol2[x_] := (cs2 . ansatz/.csNew)/.X->x;
		map1New[y_]:= map1/.x->y;
		map2New[y_]:= map2/.x->y;
		sol1New[x_]:= sol1[map1New[x]];
		sol2New[x_]:= sol2[map2New[x]];
		poly1 = Function[\[Sigma],Which[\[Sigma]<\[Sigma]p,sol1New[\[Sigma]],\[Sigma]>\[Sigma]p,sol2New[\[Sigma]]]]; 
		poly2 = Function[r, Which[r<r0,Z[2/r,\[Xi]]poly1[2/r],r>r0,Z[2/r,\[Xi]]poly1[2/r]]];
		Return[{poly1,poly2}]
]


(* ::Section:: *)
(*Overall Module*)


Options[ReggeWheelerHyperboloidal]={"GridPoints" -> 32};


ReggeWheelerHyperboloidal[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]]:=
 Module[{r0, M, w, grid,Xgrid,prec,npts, S, R, R\[Sigma], Z},
			
		(* Initial setup *)
		M = 1;
		r0 = orbit["p"];
		
		(* Initialising number of grid points *)
		npts = OptionValue["GridPoints"];
		
		(* Setting working precision *)
		prec = Precision[r0];
		
		(* Initialising Chebyshev-Gauss-Lobatto grid *)
		Xgrid = SetPrecision[Cos[(Range[0,npts]\[Pi])/npts][[2;;-2]],prec];
		
		(* Output *)
		R = HyperboloidalSolver[r0, l, m, Xgrid,
				"GridPoints" -> npts
			][[2]];
			
		R\[Sigma] = HyperboloidalSolver[r0, l, m, Xgrid,
				"GridPoints" -> npts
			][[1]];
		
		S = SpinWeightedSpheroidalHarmonicS[s, l, m, 0];
		Z = <| "\[ScriptCapitalI]" -> R\[Sigma][0], "\[ScriptCapitalH]" -> R\[Sigma][1] |>;
		w = \[Omega][\[CapitalOmega][r0,M],m];
	
		assoc = <|  "s" -> 2,
					"l" -> l,
					"m" -> m,
					"\[Omega]" -> w,
					"Type" -> {"PointParticleCircular","Orbital Radius"->ToString[r0] <>"M"},
					"RadialFunction" -> R,
					"AngularFunction" -> S,
					"Amplitudes" -> Z,
					"Method" -> {"Hyperboloidal", "GridPoints"->npts}
				|>;
			
	ReggeWheeler`ReggeWheelerMode`ReggeWheelerMode[assoc]
]


(* ::Section:: *)
(*End Package*)


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];

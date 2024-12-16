(* ::Package:: *)

(* ::Chapter:: *)
(*Hyperboloidal*)


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
(*Unprotect symbols*)


(*ClearAttributes[{ReggeWheelerHyperboloidalMode, ReggeWheelerHyperboloidal}, {Protected}];*)


(* ::Subsection:: *)
(*Usage messages*)


(*ReggeWheelerHyperboloidalMode::usage = "ReggeWheelerHyperboloidalMode[assoc] is an object which represents a Regge Wheeler mode obtained using hyperboloidal compactification.";
ReggeWheelerHyperboloidal::usage = "ReggeWheelerHyperboloidal[s, l, m, n, orbit] produces a "<>
	"ReggeWheelerHyperboloidalMode representing an inhomogeneous solution to the Regge-Wheeler equation with a point particle source constructed using hyperboloidal compactification. The formulation is based on arXiv:gr-qc/2202.01794 and arXiv:gr-qc/2411.14976."*)


(* ::Subsection:: *)
(*Error Messages*)


(*ReggeWheelerHyperboloidal::nospin = "Regge-Wheeler perturbations are only for Schwarzschild black holes but spin `1` is not zero."
ReggeWheelerHyperboloidal::eccentricity = "This package does not currently accept eccentric orbits. Please set eccentricity ('e') to zero."
ReggeWheelerHyperboloidal::spin2field = "This package currently only works for spin = 2 fields, but the fluxes and radial fns. are correct for spin = -2. Please set field spin ('s') to two."
ReggeWheelerHyperboloidal::inclination = "This package currently only works for orbits in the equatorial plane. Please set orbital inclination ('x') to one."
ReggeWheelerHyperboloidal::eccentricitymode = "This package currently only works for circular orbit modes ('m'). Please set the eccentricity mode ('n') to zero."*)


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
	
(* coefficient of height fn in exponent  *)
	\[Xi][r_,m_,M_]:= -I \[Omega][\[CapitalOmega][r,M],m] 4M;
	
(* Coefficients of hyperboloidal master eqn. operator *)
	\[Alpha]2[\[Sigma]_] := \[Sigma]^2 (1-\[Sigma]);
	\[Alpha]1[\[Sigma]_,s_]:= \[Sigma](2-3\[Sigma])+s(1-2\[Sigma]^2);
	(* Distinct coefficients for the even/odd parity equations *)
	\[Alpha]0Odd[\[Sigma]_,l_,M_,s_]:= -(s^2 (1+\[Sigma])+2s \[Sigma]+(4M^2 VeffOdd[2/\[Sigma],l,M])/((1-\[Sigma])\[Sigma]^2)); 
	\[Alpha]0Even[\[Sigma]_,l_,M_,s_]:= -(s^2 (1+\[Sigma])+2s \[Sigma]+(4M^2 VeffEven[2/\[Sigma],l,M])/((1-\[Sigma])\[Sigma]^2));


(* ::Section:: *)
(*Spectral ODE Solver*)


Options[HyperboloidalSolver]={"WorkingPrecision"->32, "GridPoints"->32};


HyperboloidalSolver[r0_, l_, m_, Xgrid_, opts:OptionsPattern[]]:=Module[
	{npts, M, \[Theta], s, \[Sigma]p, prec, map1, map2, InvMap1, InvMap2, x, X, \[Phi], S1, S2, ansatz, Dansatz, D2ansatz, Dmap1, Dmap2, 
	cs1, cs2, cs, DH, A, B, ansatz1D1, ansatz1D2, ansatz2D1, ansatz2D2, \[Alpha]21, 
	\[Alpha]22, \[Alpha]11, \[Alpha]12, \[Alpha]01, \[Alpha]02, BCs, BCsRHS, D\[Alpha]2, ODEs, juncs, juncs1,
	juncs2, fill, juncsRHS, dom1, dom2, Mat, LARHS, 
	sols, sols2, csNew, sol1, sol2, map1New, map2New, y, sol1New, sol2New, poly1, poly2},
	(*Module internal number of points on grid (maybe unnecessary)*)
		npts = OptionValue["GridPoints"];
		
	(* Initial setup *)
		M = 1;
		\[Theta] = \[Pi]/2;
		s = \[Xi][r0,m,M];
		(* Source radial position in hyp coords *)
		\[Sigma]p = 2/r0;
		(* Setting working precision *)
		prec = OptionValue["WorkingPrecision"];
		
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
		{A,B} = { (2E^(-s H[\[Sigma]p]))/(1-\[Sigma]p) (2 S1+\[Sigma]p^2/(1-\[Sigma]p) (1-s(1-\[Sigma]p)(DH))S2),-((2 \[Sigma]p^2 E^(-s H[\[Sigma]p]))/(1-\[Sigma]p))S2};
		
	(* Transforming ansatz derivatives  *)	
		{ansatz1D1, ansatz2D1, ansatz1D2, ansatz2D2} = {Dmap1*Dansatz, Dmap2*Dansatz, Dmap1^2*D2ansatz, Dmap2^2*D2ansatz};
		
	(* Obtaining master fn operator coefficients on Chebyshev-Gauss-Lobatto grid *)	
		{\[Alpha]21, \[Alpha]22, \[Alpha]11, \[Alpha]12, \[Alpha]01, \[Alpha]02} = If[EvenQ[l+m],
										{\[Alpha]2[InvMap1], \[Alpha]2[InvMap2], \[Alpha]1[InvMap1,s], \[Alpha]1[InvMap2,s], \[Alpha]0Even[InvMap1,l,M,s], \[Alpha]0Even[InvMap2,l,M,s]},
										{\[Alpha]2[InvMap1], \[Alpha]2[InvMap2], \[Alpha]1[InvMap1,s], \[Alpha]1[InvMap2,s], \[Alpha]0Odd[InvMap1,l,M,s], \[Alpha]0Odd[InvMap2,l,M,s]}];
		
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
		juncsRHS = {B/\[Alpha]2[\[Sigma]p ],A/\[Alpha]2[\[Sigma]p ]+B (D\[Alpha]2 -\[Alpha]1[\[Sigma]p,s ])/\[Alpha]2[\[Sigma]p ]^2};
		
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
		poly2 = Function[r,Which[2/r<\[Sigma]p,sol1New[2/r],2/r>\[Sigma]p,sol2New[2/r]]];
		Return[{poly1,poly2}]
]


(* ::Section:: *)
(*Overall Module*)


Options[ReggeWheelerHyperboloidal]={"WorkingPrecision" -> 32, "GridPoints" -> 32};


ReggeWheelerHyperboloidal[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]]:=
 Module[{r0, M, w, grid,Xgrid,prec,npts, S, R, R\[Sigma], Z},
	  
	  (* Error messages *)
	  If[orbit["a"] != 0,
	    Message[ReggeWheelerHyperboloidal::nospin, orbit["a"]];
	    Return[$Failed];
		];
		
	  If[orbit["e"] != 0,
	    Message[ReggeWheelerHyperboloidal::eccentricity];
	    Return[$Failed];
		];
		
	  If[s != 2,
	    Message[ReggeWheelerHyperboloidal::spin2field];
	    Return[$Failed];
		];
		
	  If[orbit["Inclination"] != 1,
	    Message[ReggeWheelerHyperboloidal::inclination];
	    Return[$Failed];
		];
		
	  If[n != 0,
	    Message[ReggeWheelerHyperboloidal::eccentricitymode];
	    Return[$Failed];
		];
			
		(* Initial setup *)
		M = 1;
		r0 = orbit["p"];
		
		(* Defining working precision *)
		prec = OptionValue["WorkingPrecision"];
		
		(* Initialising number of grid points *)
		npts = OptionValue["GridPoints"];
		
		(* Initialising Chebyshev-Gauss-Lobatto grid *)
		Xgrid = SetPrecision[Cos[(Range[0,npts]\[Pi])/npts][[2;;-2]],prec];
		
		(* Output *)
		R = HyperboloidalSolver[r0, l, m, Xgrid,
				"WorkingPrecision"->prec,
				"GridPoints" -> npts
			][[2]];
			
		R\[Sigma] = HyperboloidalSolver[r0, l, m, Xgrid,
				"WorkingPrecision"->prec,
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


(* ::Section::Closed:: *)
(*ReggeWheelerMode*)


(* ::Subsection::Closed:: *)
(*Output format*)


(*ReggeWheelerHyperboloidalMode /:
 MakeBoxes[rwm:ReggeWheelerHyperboloidalMode[assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", assoc["s"]}], "  ",
                  BoxForm`SummaryItem[{"l: ", assoc["l"]}], "  ",
                  BoxForm`SummaryItem[{"m: ", assoc["m"]}], "  ",
                  BoxForm`SummaryItem[{"\[Omega]: ", assoc["\[Omega]"]}]}],
             BoxForm`SummaryItem[{"Type: ", First[assoc["Type"]]}]};
  extended = {BoxForm`SummaryItem[{"Amplitude at \[ScriptCapitalI]: ", assoc["Amplitudes"][[1]]}],
              BoxForm`SummaryItem[{"Amplitude at \[ScriptCapitalH]: ", assoc["Amplitudes"][[2]]}],
              BoxForm`SummaryItem[{"Type details: ", Column[Rest[assoc["Type"]]]}]};
  BoxForm`ArrangeSummaryBox[
    ReggeWheelerHyperboloidalMode,
    rwm,
    None,
    summary,
    extended,
    form
  ]
];*)


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


(*ReggeWheelerHyperboloidalMode[assoc_]["EnergyFlux"] := EnergyFlux[ReggeWheelerHyperboloidalMode[assoc]];*)


(*ReggeWheelerHyperboloidalMode[assoc_]["AngularMomentumFlux"] := AngularMomentumFlux[ReggeWheelerHyperboloidalMode[assoc]];*)


(*ReggeWheelerHyperboloidalMode[assoc_]["Fluxes"] := <|"Energy" -> ReggeWheelerHyperboloidalMode[assoc]["EnergyFlux"], "AngularMomentum" -> ReggeWheelerHyperboloidalMode[assoc]["AngularMomentumFlux"]|>;*)


(*ReggeWheelerHyperboloidalMode[assoc_][string_] := assoc[string];*)


(* ::Section::Closed:: *)
(*Fluxes*)


(* ::Subsection:: *)
(*Energy Flux*)


(*EnergyFlux[mode_ReggeWheelerHyperboloidalMode] :=
	Module[{l, m, \[Xi], Z, FluxInf, FluxH, r},
	l = mode["l"];
	m = mode["m"];
	Z = mode["Amplitudes"];
	\[Xi] = -I mode["\[Omega]"] 4;
	
	
	FluxInf = ((l+2)!/(l-2)!)1/(256\[Pi]*16)Abs[\[Xi] Z[[1]]]^2;
	FluxH = ((l+2)!/(l-2)!)1/(256\[Pi]*16)Abs[\[Xi] Z[[2]]]^2;
	
	<| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxH |>
];*)


(* ::Subsection:: *)
(*Angular Momentum Flux*)


(*AngularMomentumFlux[mode_ReggeWheelerHyperboloidalMode] :=
	Module[{l, m, \[Xi], Z, FluxInf, FluxH, r},
	l = mode["l"];
	m = mode["m"];
	Z = mode["Amplitudes"];
	\[Xi] = -I mode["\[Omega]"] 4;
	
	
	FluxInf = I*m*\[Xi]((l+2)!/(l-2)!)1/(64\[Pi]*16)Abs[Z[[1]]]^2;
	FluxH = I*m*\[Xi]((l+2)!/(l-2)!)1/(64\[Pi]*16)Abs[Z[[2]]]^2;
	
	<| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxH |>
];*)


(* ::Section:: *)
(*End Package*)


(* ::Subsection:: *)
(*Protect symbols*)


(*SetAttributes[{ReggeWheelerHyperboloidalMode, ReggeWheelerHyperboloidal}, {Protected, ReadProtected}];*)


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];

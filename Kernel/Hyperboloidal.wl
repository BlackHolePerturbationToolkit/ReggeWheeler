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


(* ::Section::Closed:: *)
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

(* RWZ gauge source terms. These differ slightly from the papers cited above to ensure results agree with the conventions used in the overall ReggeWheeler package *)
	SourceTerm1Even[r_,M_,l_,m_,\[Theta]_]:=(l(l+1))/8 (( p[r]q[r,l,M] \[ScriptCapitalE][r,M])/(r f[r,M]\[CapitalLambda][r,l,M]) (L0[r,M]^2/\[ScriptCapitalE][r,M]^2 f[r,M]^2 \[CapitalLambda][r,l,M]-(\[Lambda][l](\[Lambda][l]+1)r^2+6\[Lambda][l]M r +15M^2))
									Conjugate[SphericalHarmonicY[l,m,\[Theta],0]]-(4p[r]L0[r,M]^2 f[r,M]^2)/(r \[ScriptCapitalE][r,M]) (l-2)!/(l+2)! Y\[Phi][l,m,\[Theta]]);
	SourceTerm2Even[r_,M_,l_,m_,\[Theta]_]:=(l(l+1))/8 (p[r] q[r,l,M]r^2 \[ScriptCapitalE][r,M])Conjugate[SphericalHarmonicY[l,m,\[Theta],0]];
	SourceTerm1Odd[r_,M_,l_,m_,\[Theta]_]:= ((l-1)(l+2))/4 ((2 p[r]f[r,M]L0[r,M])/(\[Lambda][l]*L[l]))X\[Theta][l,m,\[Theta]];
	SourceTerm2Odd[r_,M_,l_,m_,\[Theta]_]:= -(((l-1)(l+2))/4)(2 p[r] r f[r,M]^2 L0[r,M])/(\[Lambda][l]*L[l])X\[Theta][l,m,\[Theta]];

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
	
(* Analytic mesh refinement (AnMR) \[Kappa] and coordinate transformation *)
	\[Kappa][r0_]:=1/2 Log[r0];
	AnMRTransform[\[Chi]_,r0_,prec_]:=N[-1(1-(2Sinh[\[Kappa][r0](1+ \[Chi])])/Sinh[2\[Kappa][r0]]),prec];


(* ::Section::Closed:: *)
(*Creating Working Precision Interpolator*)


	(* Data from results obtained from our own tests on accuracy of results based on number of Chebyshev nodes and working precision at various radii up to r0 = 10^4M and modes up to lmax = 15 *)
	InterpolatorData = {{1.`,2.`,1.`,26.`},{1.`,2.`,2.`,25.`},{1.`,3.`,1.`,28.`},{1.`,3.`,2.`,28.`},{1.`,3.`,3.`,27.`},{1.`,4.`,1.`,31.`},{1.`,4.`,2.`,30.`},{1.`,4.`,3.`,30.`},{1.`,4.`,4.`,29.`},{1.`,5.`,1.`,34.`},{1.`,5.`,2.`,33.`},{1.`,5.`,3.`,33.`},{1.`,5.`,4.`,33.`},{1.`,5.`,5.`,32.`},{1.`,6.`,1.`,37.`},{1.`,6.`,2.`,36.`},{1.`,6.`,3.`,36.`},{1.`,6.`,4.`,35.`},{1.`,6.`,5.`,35.`},{1.`,6.`,6.`,34.`},{1.`,7.`,1.`,39.`},{1.`,7.`,2.`,39.`},{1.`,7.`,3.`,38.`},{1.`,7.`,4.`,38.`},{1.`,7.`,5.`,37.`},{1.`,7.`,6.`,37.`},{1.`,7.`,7.`,36.`},{1.`,8.`,1.`,43.`},{1.`,8.`,2.`,41.`},{1.`,8.`,3.`,41.`},{1.`,8.`,4.`,40.`},{1.`,8.`,5.`,40.`},{1.`,8.`,6.`,39.`},{1.`,8.`,7.`,40.`},{1.`,8.`,8.`,39.`},{1.`,9.`,1.`,47.`},{1.`,9.`,2.`,44.`},{1.`,9.`,3.`,43.`},{1.`,9.`,4.`,43.`},{1.`,9.`,5.`,42.`},{1.`,9.`,6.`,43.`},{1.`,9.`,7.`,42.`},{1.`,9.`,8.`,42.`},{1.`,9.`,9.`,41.`},{1.`,10.`,1.`,51.`},{1.`,10.`,2.`,46.`},{1.`,10.`,3.`,46.`},{1.`,10.`,4.`,45.`},{1.`,10.`,5.`,45.`},{1.`,10.`,6.`,45.`},{1.`,10.`,7.`,45.`},{1.`,10.`,8.`,44.`},{1.`,10.`,9.`,44.`},{1.`,10.`,10.`,43.`},{1.`,11.`,1.`,54.`},{1.`,11.`,2.`,49.`},{1.`,11.`,3.`,48.`},{1.`,11.`,4.`,48.`},{1.`,11.`,5.`,47.`},{1.`,11.`,6.`,48.`},{1.`,11.`,7.`,47.`},{1.`,11.`,8.`,47.`},{1.`,11.`,9.`,46.`},{1.`,11.`,10.`,47.`},{1.`,11.`,11.`,46.`},{1.`,12.`,1.`,58.`},{1.`,12.`,2.`,51.`},{1.`,12.`,3.`,51.`},{1.`,12.`,4.`,50.`},{1.`,12.`,5.`,50.`},{1.`,12.`,6.`,50.`},{1.`,12.`,7.`,50.`},{1.`,12.`,8.`,49.`},{1.`,12.`,9.`,49.`},{1.`,12.`,10.`,49.`},{1.`,12.`,11.`,49.`},{1.`,12.`,12.`,48.`},{1.`,13.`,1.`,62.`},{1.`,13.`,2.`,54.`},{1.`,13.`,3.`,53.`},{1.`,13.`,4.`,53.`},{1.`,13.`,5.`,53.`},{1.`,13.`,6.`,53.`},{1.`,13.`,7.`,52.`},{1.`,13.`,8.`,52.`},{1.`,13.`,9.`,51.`},{1.`,13.`,10.`,52.`},{1.`,13.`,11.`,51.`},{1.`,13.`,12.`,51.`},{1.`,13.`,13.`,50.`},{1.`,14.`,1.`,66.`},{1.`,14.`,2.`,57.`},{1.`,14.`,3.`,56.`},{1.`,14.`,4.`,55.`},{1.`,14.`,5.`,56.`},{1.`,14.`,6.`,55.`},{1.`,14.`,7.`,55.`},{1.`,14.`,8.`,54.`},{1.`,14.`,9.`,55.`},{1.`,14.`,10.`,54.`},{1.`,14.`,11.`,54.`},{1.`,14.`,12.`,53.`},{1.`,14.`,13.`,54.`},{1.`,14.`,14.`,52.`},{1.`,15.`,1.`,70.`},{1.`,15.`,2.`,61.`},{1.`,15.`,3.`,58.`},{1.`,15.`,4.`,58.`},{1.`,15.`,5.`,58.`},{1.`,15.`,6.`,58.`},{1.`,15.`,7.`,57.`},{1.`,15.`,8.`,57.`},{1.`,15.`,9.`,57.`},{1.`,15.`,10.`,57.`},{1.`,15.`,11.`,56.`},{1.`,15.`,12.`,56.`},{1.`,15.`,13.`,55.`},{1.`,15.`,14.`,56.`},{1.`,15.`,15.`,55.`},{2.`,2.`,1.`,36.`},{2.`,2.`,2.`,34.`},{2.`,3.`,1.`,40.`},{2.`,3.`,2.`,40.`},{2.`,3.`,3.`,39.`},{2.`,4.`,1.`,46.`},{2.`,4.`,2.`,44.`},{2.`,4.`,3.`,45.`},{2.`,4.`,4.`,43.`},{2.`,5.`,1.`,49.`},{2.`,5.`,2.`,50.`},{2.`,5.`,3.`,48.`},{2.`,5.`,4.`,50.`},{2.`,5.`,5.`,48.`},{2.`,6.`,1.`,55.`},{2.`,6.`,2.`,53.`},{2.`,6.`,3.`,54.`},{2.`,6.`,4.`,53.`},{2.`,6.`,5.`,54.`},{2.`,6.`,6.`,52.`},{2.`,7.`,1.`,59.`},{2.`,7.`,2.`,59.`},{2.`,7.`,3.`,58.`},{2.`,7.`,4.`,59.`},{2.`,7.`,5.`,57.`},{2.`,7.`,6.`,59.`},{2.`,7.`,7.`,57.`},{2.`,8.`,1.`,65.`},{2.`,8.`,2.`,63.`},{2.`,8.`,3.`,64.`},{2.`,8.`,4.`,62.`},{2.`,8.`,5.`,63.`},{2.`,8.`,6.`,62.`},{2.`,8.`,7.`,63.`},{2.`,8.`,8.`,61.`},{2.`,9.`,1.`,68.`},{2.`,9.`,2.`,69.`},{2.`,9.`,3.`,67.`},{2.`,9.`,4.`,68.`},{2.`,9.`,5.`,67.`},{2.`,9.`,6.`,68.`},{2.`,9.`,7.`,66.`},{2.`,9.`,8.`,68.`},{2.`,9.`,9.`,66.`},{2.`,10.`,1.`,74.`},{2.`,10.`,2.`,72.`},{2.`,10.`,3.`,73.`},{2.`,10.`,4.`,71.`},{2.`,10.`,5.`,73.`},{2.`,10.`,6.`,71.`},{2.`,10.`,7.`,72.`},{2.`,10.`,8.`,71.`},{2.`,10.`,9.`,72.`},{2.`,10.`,10.`,70.`},{2.`,11.`,1.`,77.`},{2.`,11.`,2.`,78.`},{2.`,11.`,3.`,76.`},{2.`,11.`,4.`,78.`},{2.`,11.`,5.`,76.`},{2.`,11.`,6.`,77.`},{2.`,11.`,7.`,76.`},{2.`,11.`,8.`,77.`},{2.`,11.`,9.`,75.`},{2.`,11.`,10.`,77.`},{2.`,11.`,11.`,75.`},{2.`,12.`,1.`,83.`},{2.`,12.`,2.`,81.`},{2.`,12.`,3.`,82.`},{2.`,12.`,4.`,81.`},{2.`,12.`,5.`,82.`},{2.`,12.`,6.`,80.`},{2.`,12.`,7.`,82.`},{2.`,12.`,8.`,80.`},{2.`,12.`,9.`,82.`},{2.`,12.`,10.`,80.`},{2.`,12.`,11.`,82.`},{2.`,12.`,12.`,79.`},{2.`,13.`,1.`,87.`},{2.`,13.`,2.`,87.`},{2.`,13.`,3.`,86.`},{2.`,13.`,4.`,87.`},{2.`,13.`,5.`,85.`},{2.`,13.`,6.`,86.`},{2.`,13.`,7.`,85.`},{2.`,13.`,8.`,86.`},{2.`,13.`,9.`,85.`},{2.`,13.`,10.`,86.`},{2.`,13.`,11.`,84.`},{2.`,13.`,12.`,86.`},{2.`,13.`,13.`,84.`},{2.`,14.`,1.`,93.`},{2.`,14.`,2.`,91.`},{2.`,14.`,3.`,92.`},{2.`,14.`,4.`,90.`},{2.`,14.`,5.`,91.`},{2.`,14.`,6.`,90.`},{2.`,14.`,7.`,91.`},{2.`,14.`,8.`,89.`},{2.`,14.`,9.`,91.`},{2.`,14.`,10.`,89.`},{2.`,14.`,11.`,91.`},{2.`,14.`,12.`,89.`},{2.`,14.`,13.`,91.`},{2.`,14.`,14.`,88.`},{2.`,15.`,1.`,96.`},{2.`,15.`,2.`,97.`},{2.`,15.`,3.`,95.`},{2.`,15.`,4.`,96.`},{2.`,15.`,5.`,94.`},{2.`,15.`,6.`,96.`},{2.`,15.`,7.`,94.`},{2.`,15.`,8.`,95.`},{2.`,15.`,9.`,94.`},{2.`,15.`,10.`,95.`},{2.`,15.`,11.`,94.`},{2.`,15.`,12.`,95.`},{2.`,15.`,13.`,93.`},{2.`,15.`,14.`,95.`},{2.`,15.`,15.`,93.`},{3.`,2.`,1.`,46.`},{3.`,2.`,2.`,43.`},{3.`,3.`,1.`,51.`},{3.`,3.`,2.`,52.`},{3.`,3.`,3.`,50.`},{3.`,4.`,1.`,60.`},{3.`,4.`,2.`,57.`},{3.`,4.`,3.`,59.`},{3.`,4.`,4.`,56.`},{3.`,5.`,1.`,64.`},{3.`,5.`,2.`,66.`},{3.`,5.`,3.`,63.`},{3.`,5.`,4.`,66.`},{3.`,5.`,5.`,63.`},{3.`,6.`,1.`,73.`},{3.`,6.`,2.`,70.`},{3.`,6.`,3.`,72.`},{3.`,6.`,4.`,70.`},{3.`,6.`,5.`,72.`},{3.`,6.`,6.`,69.`},{3.`,7.`,1.`,78.`},{3.`,7.`,2.`,80.`},{3.`,7.`,3.`,77.`},{3.`,7.`,4.`,79.`},{3.`,7.`,5.`,76.`},{3.`,7.`,6.`,79.`},{3.`,7.`,7.`,76.`},{3.`,8.`,1.`,87.`},{3.`,8.`,2.`,84.`},{3.`,8.`,3.`,86.`},{3.`,8.`,4.`,83.`},{3.`,8.`,5.`,85.`},{3.`,8.`,6.`,83.`},{3.`,8.`,7.`,85.`},{3.`,8.`,8.`,82.`},{3.`,9.`,1.`,91.`},{3.`,9.`,2.`,93.`},{3.`,9.`,3.`,90.`},{3.`,9.`,4.`,92.`},{3.`,9.`,5.`,90.`},{3.`,9.`,6.`,92.`},{3.`,9.`,7.`,89.`},{3.`,9.`,8.`,92.`},{3.`,9.`,9.`,89.`},{3.`,10.`,1.`,100.`},{3.`,10.`,2.`,97.`},{3.`,10.`,3.`,99.`},{3.`,10.`,4.`,97.`},{3.`,10.`,5.`,99.`},{3.`,10.`,6.`,96.`},{3.`,10.`,7.`,99.`},{3.`,10.`,8.`,96.`},{3.`,10.`,9.`,98.`},{3.`,10.`,10.`,95.`},{3.`,11.`,1.`,104.`},{3.`,11.`,2.`,106.`},{3.`,11.`,3.`,103.`},{3.`,11.`,4.`,106.`},{3.`,11.`,5.`,103.`},{3.`,11.`,6.`,105.`},{3.`,11.`,7.`,103.`},{3.`,11.`,8.`,105.`},{3.`,11.`,9.`,102.`},{3.`,11.`,10.`,105.`},{3.`,11.`,11.`,102.`},{3.`,12.`,1.`,113.`},{3.`,12.`,2.`,110.`},{3.`,12.`,3.`,112.`},{3.`,12.`,4.`,110.`},{3.`,12.`,5.`,112.`},{3.`,12.`,6.`,109.`},{3.`,12.`,7.`,112.`},{3.`,12.`,8.`,109.`},{3.`,12.`,9.`,112.`},{3.`,12.`,10.`,109.`},{3.`,12.`,11.`,112.`},{3.`,12.`,12.`,108.`},{3.`,13.`,1.`,118.`},{3.`,13.`,2.`,119.`},{3.`,13.`,3.`,117.`},{3.`,13.`,4.`,119.`},{3.`,13.`,5.`,116.`},{3.`,13.`,6.`,119.`},{3.`,13.`,7.`,116.`},{3.`,13.`,8.`,118.`},{3.`,13.`,9.`,116.`},{3.`,13.`,10.`,118.`},{3.`,13.`,11.`,115.`},{3.`,13.`,12.`,118.`},{3.`,13.`,13.`,115.`},{3.`,14.`,1.`,127.`},{3.`,14.`,2.`,124.`},{3.`,14.`,3.`,126.`},{3.`,14.`,4.`,123.`},{3.`,14.`,5.`,125.`},{3.`,14.`,6.`,123.`},{3.`,14.`,7.`,125.`},{3.`,14.`,8.`,122.`},{3.`,14.`,9.`,125.`},{3.`,14.`,10.`,122.`},{3.`,14.`,11.`,125.`},{3.`,14.`,12.`,122.`},{3.`,14.`,13.`,125.`},{3.`,14.`,14.`,122.`},{3.`,15.`,1.`,131.`},{3.`,15.`,2.`,133.`},{3.`,15.`,3.`,130.`},{3.`,15.`,4.`,132.`},{3.`,15.`,5.`,130.`},{3.`,15.`,6.`,132.`},{3.`,15.`,7.`,129.`},{3.`,15.`,8.`,132.`},{3.`,15.`,9.`,129.`},{3.`,15.`,10.`,131.`},{3.`,15.`,11.`,129.`},{3.`,15.`,12.`,131.`},{3.`,15.`,13.`,129.`},{3.`,15.`,14.`,131.`},{3.`,15.`,15.`,128.`},{4.`,2.`,1.`,56.`},{4.`,2.`,2.`,52.`},{4.`,3.`,1.`,62.`},{4.`,3.`,2.`,64.`},{4.`,3.`,3.`,61.`},{4.`,4.`,1.`,74.`},{4.`,4.`,2.`,70.`},{4.`,4.`,3.`,73.`},{4.`,4.`,4.`,69.`},{4.`,5.`,1.`,79.`},{4.`,5.`,2.`,82.`},{4.`,5.`,3.`,78.`},{4.`,5.`,4.`,82.`},{4.`,5.`,5.`,78.`},{4.`,6.`,1.`,91.`},{4.`,6.`,2.`,87.`},{4.`,6.`,3.`,90.`},{4.`,6.`,4.`,87.`},{4.`,6.`,5.`,90.`},{4.`,6.`,6.`,86.`},{4.`,7.`,1.`,97.`},{4.`,7.`,2.`,100.`},{4.`,7.`,3.`,96.`},{4.`,7.`,4.`,99.`},{4.`,7.`,5.`,95.`},{4.`,7.`,6.`,99.`},{4.`,7.`,7.`,95.`},{4.`,8.`,1.`,109.`},{4.`,8.`,2.`,105.`},{4.`,8.`,3.`,108.`},{4.`,8.`,4.`,104.`},{4.`,8.`,5.`,107.`},{4.`,8.`,6.`,104.`},{4.`,8.`,7.`,107.`},{4.`,8.`,8.`,103.`},{4.`,9.`,1.`,114.`},{4.`,9.`,2.`,117.`},{4.`,9.`,3.`,113.`},{4.`,9.`,4.`,116.`},{4.`,9.`,5.`,113.`},{4.`,9.`,6.`,116.`},{4.`,9.`,7.`,112.`},{4.`,9.`,8.`,116.`},{4.`,9.`,9.`,112.`},{4.`,10.`,1.`,126.`},{4.`,10.`,2.`,122.`},{4.`,10.`,3.`,125.`},{4.`,10.`,4.`,122.`},{4.`,10.`,5.`,125.`},{4.`,10.`,6.`,121.`},{4.`,10.`,7.`,125.`},{4.`,10.`,8.`,121.`},{4.`,10.`,9.`,125.`},{4.`,10.`,10.`,120.`},{4.`,11.`,1.`,131.`},{4.`,11.`,2.`,134.`},{4.`,11.`,3.`,130.`},{4.`,11.`,4.`,134.`},{4.`,11.`,5.`,130.`},{4.`,11.`,6.`,133.`},{4.`,11.`,7.`,130.`},{4.`,11.`,8.`,133.`},{4.`,11.`,9.`,129.`},{4.`,11.`,10.`,133.`},{4.`,11.`,11.`,129.`},{4.`,12.`,1.`,143.`},{4.`,12.`,2.`,139.`},{4.`,12.`,3.`,143.`},{4.`,12.`,4.`,139.`},{4.`,12.`,5.`,142.`},{4.`,12.`,6.`,138.`},{4.`,12.`,7.`,142.`},{4.`,12.`,8.`,138.`},{4.`,12.`,9.`,142.`},{4.`,12.`,10.`,138.`},{4.`,12.`,11.`,142.`},{4.`,12.`,12.`,137.`},{4.`,13.`,1.`,149.`},{4.`,13.`,2.`,152.`},{4.`,13.`,3.`,148.`},{4.`,13.`,4.`,151.`},{4.`,13.`,5.`,147.`},{4.`,13.`,6.`,151.`},{4.`,13.`,7.`,147.`},{4.`,13.`,8.`,150.`},{4.`,13.`,9.`,147.`},{4.`,13.`,10.`,150.`},{4.`,13.`,11.`,146.`},{4.`,13.`,12.`,150.`},{4.`,13.`,13.`,146.`},{4.`,14.`,1.`,161.`},{4.`,14.`,2.`,157.`},{4.`,14.`,3.`,160.`},{4.`,14.`,4.`,156.`},{4.`,14.`,5.`,159.`},{4.`,14.`,6.`,156.`},{4.`,14.`,7.`,159.`},{4.`,14.`,8.`,155.`},{4.`,14.`,9.`,159.`},{4.`,14.`,10.`,155.`},{4.`,14.`,11.`,159.`},{4.`,14.`,12.`,155.`},{4.`,14.`,13.`,159.`},{4.`,14.`,14.`,155.`},{4.`,15.`,1.`,166.`},{4.`,15.`,2.`,169.`},{4.`,15.`,3.`,165.`},{4.`,15.`,4.`,168.`},{4.`,15.`,5.`,165.`},{4.`,15.`,6.`,168.`},{4.`,15.`,7.`,164.`},{4.`,15.`,8.`,168.`},{4.`,15.`,9.`,164.`},{4.`,15.`,10.`,167.`},{4.`,15.`,11.`,164.`},{4.`,15.`,12.`,167.`},{4.`,15.`,13.`,164.`},{4.`,15.`,14.`,167.`},{4.`,15.`,15.`,163.`}};
	
	(* Building an interpolating function from above data *)
	precFit = Interpolation[Table[{InterpolatorData[[i,1;;3]],InterpolatorData[[i,4]]},{i,Length[InterpolatorData]}], InterpolationOrder->1, Method->"Hermite"];
	
	(* Stopping extrapolation by mapping all variables to the unit interval and placing it on the unit sphere *)
	maxRangeInterp={4,15,15};
	testPoint={3.2,20,1};
	cluster=Transpose[Transpose[InterpolatorData][[1;;3]]];
	rISCO=6;
	map1[r_]=a r +b/.Solve[{0==a Log10[rISCO] +b,1==a 4 +b},{a,b}][[1]];
	map2[l_]=a l +b/.Solve[{0==a 2 +b,1==a 15 +b},{a,b}][[1]];
	map3[m_]=a m +b/.Solve[{0==a 1 +b,1==a 15 +b},{a,b}][[1]];
	maps[r_,l_,m_]={map1[r],map2[l],map3[m]};
	mappedTestPointPosition=maps[testPoint[[1]],testPoint[[2]],testPoint[[3]]];
	
	(* Building extrapolator *)
	extrapolator[r_,l_,m_]:=Module[
		{input,nP,mapNP,mapInput,npPrec,newPrec},
		input={Log10[r],l,m};
		mapInput=maps[input[[1]],input[[2]],input[[3]]];
		Round[Mean[Table[
					nP=Nearest[cluster,input,3][[index]];
					mapNP=maps[nP[[1]],nP[[2]],nP[[3]]];
					npPrec=precFit[nP[[1]],nP[[2]],nP[[3]]];
					newPrec=Round[npPrec Norm[mapInput]/Norm[mapNP]]
		,{index,3}]]]
	];
	
	(* Defining final function to determine minimum necessary working precision *)
	necessaryMinPrecision[r_,l_,m_]:=If[Log10[rISCO]<= Log10[r]<= 4 && 2<=l<= 15 && 1<=m<= 15, precFit[Log10[r],l,m], extrapolator[r,l,m]]


(* ::Section::Closed:: *)
(*Coordinate Transformation*)


DomainMapping[r0_,x_,X_,prec_]:= 
		Module[{\[Sigma]p,\[CapitalSigma], \[Sigma]grid1, \[Sigma]grid2,\[Sigma]grid3, M, AB1, AB2,AB3, map1, map2,map3, InvMap1, InvMap2,InvMap3, \[Chi]ofX, \[Chi]of\[Sigma], \[Sigma]of\[Chi], a, b},
			
			(* Initial setup *)
			M = 1;
			\[Sigma]p =(2M)/r0;
			\[CapitalSigma] = If[r0< 10, 1/2, 1/MantissaExponent[r0][[2]]];
			
			\[Sigma]grid1 = {0,\[Sigma]p};
			\[Sigma]grid2 = {\[Sigma]p,\[CapitalSigma]};
			\[Sigma]grid3 = {\[CapitalSigma],1};
			
			(* Solving for the transform coefficients *);
			AB1 = Solve[{a*\[Sigma]grid1[[1]]+b==-1,a*\[Sigma]grid1[[2]]+b==1}][[1]];
			AB2 = Solve[{a*\[Sigma]grid2[[1]]+b==-1,a*\[Sigma]grid2[[2]]+b==1}][[1]];
			AB3 =Solve[{a*\[Sigma]grid3[[1]]+b==-1,a*\[Sigma]grid3[[2]]+b==1}][[1]];
			
			(* Map to hyp. coords. *);
			map1 = a*x+b/.AB1;
			map2 = a*x+b/.AB2;
			map3 = a*x+b/.AB3;
			
			(* Map to Sch. coords *);
			InvMap1 = (X-b)/a/.AB1;
			InvMap2 = (X-b)/a/.AB2;
			InvMap3 = (X-b)/a/.AB3;
			
			\[Chi]ofX = \[Chi]/.Solve[AnMRTransform[\[Chi],r0,prec]==X,\[Chi]][[1]]/.C[1]->0//Quiet;
			\[Chi]of\[Sigma] = \[Chi]ofX/.X->map2;
			\[Sigma]of\[Chi] = InvMap3/.X->AnMRTransform[X,r0,prec];
			
			Return[{\[CapitalSigma],map1, map2,map3, InvMap1, InvMap2,InvMap3, \[Chi]ofX, \[Chi]of\[Sigma], \[Sigma]of\[Chi]}]
];


(* ::Section::Closed:: *)
(*Spectral ODE Solver*)


Options[HyperboloidalSolver]={"GridPoints"->32};


HyperboloidalSolver[r0_, l_, m_, Xgrid_, opts:OptionsPattern[]]:=Module[
	{npts, M, \[Theta], \[Xi],\[CapitalSigma], \[Sigma]p, prec, DM, DM2,map1, map2,map3, InvMap1, InvMap2,InvMap3, AnMRMap\[Chi]X, AnMRMap\[Chi]\[Sigma], AnMRMap\[Sigma]\[Chi], d\[Chi]dX, d2\[Chi]dX2, x, X, \[Phi], S1, S2, ansatz, Dansatz, D2ansatz, Dmap1, Dmap2, Dmap3,
	cs1, cs2, cs3,cs, DH, A, B, ansatz1D1, ansatz1D2, ansatz2D1, ansatz2D2,ansatz3D1, ansatz3D2, \[Alpha]21, 
	\[Alpha]22,\[Alpha]23, \[Alpha]11, \[Alpha]12,\[Alpha]13, \[Alpha]01, \[Alpha]02,\[Alpha]03, BCs, BCsRHS, D\[Alpha]2, ODEs, juncs, juncs1,
	juncs2, fill, juncsRHS, dom1, dom2, dom3,Mat, LARHS, 
	sols, sols2, csNew, ansatzUneval, sol1, sol2,sol3, map1New, map2New, map3New,y, sol1New, sol2New,sol3New, poly1, poly2},
		
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
		If[prec < necessaryMinPrecision[r0,l,m], 
			prec = necessaryMinPrecision[r0,l,m];
		];
		
		(* Obtaining differentiation matrices *)
		 DM =-NDSolve`FiniteDifferenceDerivative[Derivative[1],N[Reverse[Xgrid],prec],DifferenceOrder->"Pseudospectral",PeriodicInterpolation->False]["DifferentiationMatrix"];
		 DM2 =NDSolve`FiniteDifferenceDerivative[Derivative[2],N[Reverse[Xgrid],prec],DifferenceOrder->"Pseudospectral",PeriodicInterpolation->False]["DifferentiationMatrix"];
		
		(* Coordinate mappings based off source position *)
		{\[CapitalSigma],map1[x_],map2[x_],map3[x_],InvMap1[X_],InvMap2[X_],InvMap3[X_],AnMRMap\[Chi]X[X_], AnMRMap\[Chi]\[Sigma][x_], AnMRMap\[Sigma]\[Chi][X_]} = DomainMapping[r0,x,X,prec];
		
		(* Source terms *)
		{S1,S2} = If[EvenQ[l+m],
			{SourceTerm1Even[r0,M,l,m,\[Theta]],SourceTerm2Even[r0,M,l,m,\[Theta]]},
			{SourceTerm1Odd[r0,M,l,m,\[Theta]],SourceTerm2Odd[r0,M,l,m,\[Theta]]}];
		
		(* Chebyshev polynomial as ansatz, using number of points on grid to determine length *)
		ansatz = Table[N[ChebyshevT[i,xx],prec],{i,0,Length[Xgrid]-1},{xx,Xgrid}]//Quiet;
		
		(* Derivatives of ansatz and mappings *)
		{Dansatz,D2ansatz}={Transpose[DM . ansatz],Transpose[DM2 . ansatz]} ;
		{Dmap1,Dmap2,Dmap3}={D[map1[x],x],D[map2[x],x],D[map3[x],x]};
		{d\[Chi]dX[X_], d2\[Chi]dX2[X_]} = {D[AnMRMap\[Chi]X[X],X], D[AnMRMap\[Chi]X[X],{X,2}]};
		
		(* Initialising table of weight coefficients *)
		{cs1,cs2,cs3} = {Table[Subscript[c, i, 1],{i,0,Length[Xgrid]-1}],Table[Subscript[c, i, 2],{i,0,Length[Xgrid]-1}],Table[Subscript[c, i, 3],{i,0,Length[Xgrid]-1}]};
		cs = Join[cs1,cs2,cs3];
		
		(*Derivative of height fn for convenience*)
		DH =(D[H[x],x])/.x->\[Sigma]p;
		
		(* Hyperboloidal junction condition coefficients *)
		{A,B} = { (2E^(-\[Xi] H[\[Sigma]p]))/(1-\[Sigma]p) (2 S1+\[Sigma]p^2/(1-\[Sigma]p) (1-\[Xi](1-\[Sigma]p)(DH))S2),-((2 \[Sigma]p^2 E^(-\[Xi] H[\[Sigma]p]))/(1-\[Sigma]p))S2};
		
		(* Transforming ansatz derivatives  *)	
		{ansatz1D1, ansatz2D1,ansatz3D1, ansatz1D2, ansatz2D2,ansatz3D2} = {Dmap1*Dansatz, Dmap2*Table[(d\[Chi]dX[AnMRTransform[Xgrid,r0,prec]])*Dansatz[[i]],{i,Length[Xgrid]}],Dmap3*Dansatz, Dmap1^2*D2ansatz, Dmap2^2*(Table[(d2\[Chi]dX2[AnMRTransform[Xgrid,r0,prec]])*Dansatz[[i]],{i,Length[Xgrid]}]+Table[(d\[Chi]dX[AnMRTransform[Xgrid,r0,prec]]^2)*D2ansatz[[i]],{i,Length[Xgrid]}]),Dmap3^2*D2ansatz};
		
		(* Obtaining master fn operator coefficients on Chebyshev-Gauss-Lobatto grid *)	
		{\[Alpha]21[X_],\[Alpha]22[X_],\[Alpha]23[X_],\[Alpha]11[X_],\[Alpha]12[X_],\[Alpha]13[X_],\[Alpha]01[X_],\[Alpha]02[X_],\[Alpha]03[X_]} = If[EvenQ[l+m],
						{\[Alpha]2[InvMap1[X]],\[Alpha]2[InvMap2[AnMRTransform[X,r0,prec]]],\[Alpha]2[InvMap3[X]],\[Alpha]1[InvMap1[X],\[Xi]],\[Alpha]1[InvMap2[AnMRTransform[X,r0,prec]],\[Xi]],\[Alpha]1[InvMap3[X],\[Xi]],\[Alpha]0Even[InvMap1[X],l,M,\[Xi]],\[Alpha]0Even[InvMap2[AnMRTransform[X,r0,prec]],l,M,\[Xi]],\[Alpha]0Even[InvMap3[X],l,M,\[Xi]]},
						{\[Alpha]2[InvMap1[X]],\[Alpha]2[InvMap2[AnMRTransform[X,r0,prec]]],\[Alpha]2[InvMap3[X]],\[Alpha]1[InvMap1[X],\[Xi]],\[Alpha]1[InvMap2[AnMRTransform[X,r0,prec]],\[Xi]],\[Alpha]1[InvMap3[X],\[Xi]],\[Alpha]0Odd[InvMap1[X],l,M,\[Xi]],\[Alpha]0Odd[InvMap2[AnMRTransform[X,r0,prec]],l,M,\[Xi]],\[Alpha]0Odd[InvMap3[X],l,M,\[Xi]]}];
		
		(* Defining boundary conditions, junction conditions and the ODE at every other grid point *)
		BCs={\[Alpha]11[-1]ansatz1D1[[;;,-1]]+\[Alpha]01[-1]ansatz[[;;,-1]],\[Alpha]13[1]ansatz3D1[[;;,1]]+\[Alpha]03[1] ansatz[[;;,1]]};
		ODEs={Table[\[Alpha]21[Xgrid][[2;;-2]] ansatz1D2[[i,2;;-2]]+\[Alpha]11[Xgrid][[2;;-2]] ansatz1D1[[i,2;;-2]]+\[Alpha]01[Xgrid][[2;;-2]] ansatz[[i,2;;-2]],{i,Length[ansatz1D2]}],
		Table[\[Alpha]22[Xgrid][[2;;-2]] ansatz2D2[[i,2;;-2]]+\[Alpha]12[Xgrid][[2;;-2]] ansatz2D1[[i,2;;-2]]+\[Alpha]02[Xgrid][[2;;-2]] ansatz[[i,2;;-2]],{i,Length[Xgrid]}],
		
		Table[\[Alpha]23[Xgrid][[2;;-2]] ansatz3D2[[i,2;;-2]]+\[Alpha]13[Xgrid][[2;;-2]] ansatz3D1[[i,2;;-2]]+\[Alpha]03[Xgrid][[2;;-2]] ansatz[[i,2;;-2]],{i,Length[Xgrid]}]};
		
		{BCs, ODEs, juncs1} = {BCs, ODEs, {((cs2 . ansatz[[;;,-1]])-(cs1 . ansatz[[;;,1]])),((cs2 . ansatz2D1[[;;,-1]])-(cs1 . ansatz1D1[[;;,1]])),((cs3 . ansatz[[;;,-1]])-(cs2 . ansatz[[;;,1]])),((cs3 . ansatz3D1[[;;,-1]])-(cs2 . ansatz2D1[[;;,1]]))}};
				
		 (* Obtaining values for junction conditions *)
		juncs2 = Table[Coefficient[juncs1[[j]],cs[[i]]], {j, Length[juncs1]}, {i, Length[cs]}];
		
		(* Placeholder column of zeroes for matrix construction (replaced by junction condition vectors in full matrix ) *)
		fill = ConstantArray[0, Length[ansatz]];
		
		(* Constructing matrices to be inverted for each domain *);
		dom1 = MapThread[Append, {MapThread[Prepend,{ODEs[[1]],BCs[[1]]}],fill}];
		dom2 = MapThread[Append,{MapThread[Prepend,{ODEs[[2]],fill}],fill}];
		dom3 = MapThread[Append,{MapThread[Prepend,{ODEs[[3]],fill}],BCs[[2]]}];
		
		(* Constructing overall block diagonal matrix to be inverted *)
		Mat =  ArrayFlatten[{{dom1,0,0},{0,dom2,0},{0,0,dom3}}];
		(*Setting junction conditions*)
		{Mat[[;;,Length[dom1]]], Mat[[;;,Length[dom1]+1]],Mat[[;;,2Length[dom1]]],Mat[[;;,2Length[dom1]+1]]} = juncs2;
		
		(*Initialising values for RHS of OVERALL equation Mx = b.*);
		BCsRHS = {0,0};
		D\[Alpha]2 = (D[\[Alpha]2[x],x])/.x->\[Sigma]p;
		juncsRHS = {B/\[Alpha]2[\[Sigma]p ],(A/\[Alpha]2[\[Sigma]p ]+B (D\[Alpha]2 -\[Alpha]1[\[Sigma]p,\[Xi] ])/\[Alpha]2[\[Sigma]p ]^2)};
		
		(* Constructing RHS vector *);
		LARHS = Join[{BCsRHS[[1]]},ConstantArray[0,Length[Xgrid]-2],juncsRHS,ConstantArray[0,Length[Xgrid]-2],{0},{0},ConstantArray[0,Length[Xgrid]-2],{BCsRHS[[2]]}];
		
		(*Inverting matrix to obtain weight coefficients *);
		sols = LinearSolve[Transpose[Mat],LARHS];
		
		
		(* Creating table of replacement rules for weight coefficients *);
		csNew = Table[cs[[i]]->sols[[i]],{i,1,Length[cs]}];
		
		(* Creating final polynomial *);
		ansatzUneval[x_]=Table[ChebyshevT[i,x],{i,0,Length[Xgrid]-1}];
		sol1[x_]:= cs1 . ansatzUneval[x]/.csNew;
		sol2[x_] := cs2 . ansatzUneval[x]/.csNew;
		    sol3[x_] := cs3 . ansatzUneval[x]/.csNew;
		map1New[y_]:= map1[y];
		    map2New[y_]:= AnMRMap\[Chi]\[Sigma][y];
		map3New[y_]:= map3[y];
		sol1New[x_]:= sol1[map1New[x]];
		    sol2New[x_]:= sol2[map2New[x]];
		sol3New[x_]:= sol3[map3New[x]];
		
		poly1 =  Function[\[Sigma],Which[\[Sigma]<\[Sigma]p,sol1New[\[Sigma]],\[Sigma]p<=\[Sigma]<=\[CapitalSigma],sol2New[\[Sigma]],\[Sigma]>\[CapitalSigma],sol3New[\[Sigma]]]];
		poly2 = Function[r, Z[2/r0,\[Xi]]poly1[2/r]];
		Return[{poly1,poly2}]
]


(* ::Section::Closed:: *)
(*Overall Module*)


Options[ReggeWheelerHyperboloidal]={"GridPoints" -> 32};


ReggeWheelerHyperboloidal[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]]:=
 Module[{r0, M, w, grid,Xgrid,prec,npts, S, solution, R, R\[Sigma], Z},
			
		(* Initial setup *)
		M = 1;
		r0 = orbit["p"];
		
		(* Initialising number of grid points *)
		npts = SetPrecision[Round[(Log[10, r0]+((l+Abs[l-m])/2)^(1/2))]*32,Precision[r0]];
		
		(* Checking working precision *)
		prec = Precision[r0];
		
		(* Initialising Chebyshev-Gauss-Lobatto grid *)
		Xgrid = Cos[(Range[0,npts]\[Pi])/npts];
		
		(* Output *)
		
		solution = HyperboloidalSolver[r0, l, m, Xgrid,
				"GridPoints" -> npts
			];
		
		R = solution[[2]];
			
		R\[Sigma] = solution[[1]];
		
		S = SpinWeightedSpheroidalHarmonicS[s, l, m, 0];
		(* \[ExponentialE]^(\[ImaginaryI]*4*M*w) factor is to ensure agreement with the overall ReggeWheeler package.*)
		w = \[Omega][\[CapitalOmega][r0,M],m];
		Z = <| "\[ScriptCapitalI]" -> R\[Sigma][0], "\[ScriptCapitalH]" -> E^(I*4*M*w)*R\[Sigma][1] |>;

	
		assoc = <|  "s" -> 2,
					"l" -> l,
					"m" -> m,
					"\[Omega]" -> w,
					"Eigenvalue" -> S["Eigenvalue"],
					"Type" -> {"PointParticleCircular","Orbital Radius"->ToString[r0] <>"M"},
					"RadialFunctions" -> R,
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

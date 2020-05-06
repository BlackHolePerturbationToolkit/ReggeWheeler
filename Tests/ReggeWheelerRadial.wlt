(* Mathematica Test File *)

{\[Psi]In, \[Psi]Up} = Values[ReggeWheelerRadial[2, 2, 0.1]];

(****************************************************************)
(* TeukolskyRadial                                              *)
(****************************************************************)
VerificationTest[
    \[Psi]In
    ,
    ReggeWheelerRadialFunction[2, 2, 0.1, <|
      "s" -> 2, "l" -> 2, "\[Omega]" -> 0.1, "Eigenvalue" -> \[Lambda]_, 
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> \[Nu]_},
      "BoundaryConditions" -> "In", 
      "Amplitudes" -> <|"Transmission" -> _|>,
      "Domain" -> {2, Infinity},
      "RadialFunction" ->
        ReggeWheeler`MST`MST`Private`MSTRadialIn[2, 2, 0, 0, 0.2, \[Nu]_, \[Lambda]_, _,
          {MachinePrecision, MachinePrecision/2, MachinePrecision/2}] |>
    ]
    ,
    TestID->"ReggeWheelerRadial",
    SameTest -> MatchQ
]

(****************************************************************)
(* InvalidKey                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]In["NotAKey"]
    ,
    Missing["KeyAbsent", "NotAKey"]
    ,
    TestID->"InvalidKey"
]


(****************************************************************)
(* BoundaryConditions                                           *)
(****************************************************************)
VerificationTest[
    \[Psi]In["BoundaryConditions"]
    ,
    "In"
    ,
    TestID->"BoundaryConditions"
]


(****************************************************************)
(* Method                                                       *)
(****************************************************************)
VerificationTest[
    \[Psi]In["Method"]
    ,
    {"MST", "RenormalizedAngularMomentum" -> _}
    ,
    TestID->"Method",
    SameTest -> MatchQ
]


(****************************************************************)
(* Eigenvalue                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]In["Eigenvalue"]
    ,
    0
    ,
    TestID->"Eigenvalue"
]


(****************************************************************)
(* SolutionFunctions                                            *)
(****************************************************************)
VerificationTest[
    \[Psi]In["RadialFunction"]
    ,
    ReggeWheelerRadialFunction[__]["RadialFunction"]
    ,
    TestID->"SolutionFunctions",
    SameTest -> MatchQ
]


(****************************************************************)
(* Numerical Evaluation                                         *)
(****************************************************************)
VerificationTest[
    \[Psi]In[10.0]
    ,
    93.12096963408008 + 17.576391050258476*I
    ,
    TestID->"Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Derivative Numerical Evaluation                              *)
(****************************************************************)
VerificationTest[
    \[Psi]In'[10.0]
    ,
    25.692572303913664 + 4.8480777278799865*I
    ,
    TestID->"Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Higher Derivative Numerical Evaluation                       *)
(****************************************************************)
VerificationTest[
    \[Psi]In''''[10.0]
    ,
    -0.07496035958290637 - 0.014105323216607984*I
    ,
    TestID->"Higher Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


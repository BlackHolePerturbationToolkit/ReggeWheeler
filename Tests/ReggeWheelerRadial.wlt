(* Mathematica Test File *)

\[Psi] = ReggeWheelerRadial[2, 2, 0.1];
\[Psi]In = \[Psi]["In"];

(****************************************************************)
(* TeukolskyRadial                                              *)
(****************************************************************)
VerificationTest[
    \[Psi]
    ,
    ReggeWheelerRadialFunction[2, 2, 0.1, 
     <|"Method" -> {"NumericalIntegration", "rmin" -> 4, "rmax" -> 20}, 
      "BoundaryConditions" -> {"In", "Up"}, 
      "SolutionFunctions" -> {ReggeWheeler`NumericalIntegration`Private`PsiIn[2, 2, 0.1, 4, 20], ReggeWheeler`NumericalIntegration`Private`PsiUp[2, 2, 0.1, 4, 20]}|>]
    ,
    TestID->"ReggeWheelerRadial",
    SameTest -> withinRoundoff
]

(****************************************************************)
(* InvalidKey                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]["NotAKey"]
    ,
    Missing["KeyAbsent", "NotAKey"]
    ,
    TestID->"InvalidKey"
]


(****************************************************************)
(* BoundaryConditions                                           *)
(****************************************************************)
VerificationTest[
    \[Psi]["BoundaryConditions"]
    ,
    {"In", "Up"}
    ,
    TestID->"BoundaryConditions"
]


(****************************************************************)
(* Method                                                       *)
(****************************************************************)
VerificationTest[
    \[Psi]["Method"]
    ,
    {"NumericalIntegration", "rmin" -> 4, "rmax" -> 20}
    ,
    TestID->"Method"
]


(****************************************************************)
(* SolutionFunctions                                            *)
(****************************************************************)
VerificationTest[
    \[Psi]["SolutionFunctions"]
    ,
    {ReggeWheeler`NumericalIntegration`Private`PsiIn[2, 2, 0.1, 4, 20], ReggeWheeler`NumericalIntegration`Private`PsiUp[2, 2, 0.1, 4, 20]}
    ,
    TestID->"SolutionFunctions"
]


(****************************************************************)
(* Numerical Evaluation                                         *)
(****************************************************************)
VerificationTest[
    \[Psi][10.0]
    ,
    <|"In" -> 0.0, 
     "Up" -> 0.0 |>
    ,
    TestID->"Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Derivative Numerical Evaluation                              *)
(****************************************************************)
VerificationTest[
    \[Psi]'[10.0]
    ,
    <|"In" -> 0.0, 
     "Up" -> 0.0 |>
    ,
    TestID->"Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Higher Derivative Numerical Evaluation                       *)
(****************************************************************)
VerificationTest[
    \[Psi]''''[10.0]
    ,
    <|"In" -> 0.0, 
     "Up" -> 0.0 |>
    ,
    TestID->"Higher Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Subcase                                                      *)
(****************************************************************)
VerificationTest[
    \[Psi]In
    ,
    ReggeWheelerRadialFunction[2, 2, 0.1, 
     <|"Method" -> {"NumericalIntegration", "rmin" -> 4, "rmax" -> 20}, 
      "BoundaryConditions" -> "In",
      "SolutionFunctions" -> ReggeWheeler`NumericalIntegration`Private`PsiIn[2, 2, 0.1, 4, 20]
      |>]
    ,
    TestID->"Subcase"
]

(****************************************************************)
(* Single Subcase Boundary Conditions                           *)
(****************************************************************)
VerificationTest[
    \[Psi]In["BoundaryConditions"]
    ,
    "In"
    ,
    TestID->"Single Subcase Boundary Conditions"
]

(****************************************************************)
(* Single Subcase Numerical Evaluation                          *)
(****************************************************************)
VerificationTest[
    \[Psi]In[10.0]
    ,
    0.0
    ,
    TestID->"Single Subcase Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Single Subcase Derivative Numerical Evaluation               *)
(****************************************************************)
VerificationTest[
    \[Psi]In'[10.0]
    ,
    0.0
    ,
    TestID->"Single Subcase Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]



(* ::Package:: *)

Paclet[
  Name -> "ReggeWheeler",
  Version -> "0.1.0",
  MathematicaVersion -> "10+",
  Creator -> "Black Hole Perturbation Toolkit",
  Description -> "A set of functions for computing solutions to the Regge Wheeler equation.",
  Extensions ->
  {
    { "Kernel",
      "Context" -> {
        "ReggeWheeler`",
        "ReggeWheeler`NumericalIntegration`",
        "ReggeWheeler`MST`MST`",
        "ReggeWheeler`MST`RenormalizedAngularMomentum`"
      }
    },

    { "Documentation",
      Language -> "English", 
      MainPage -> "Guides/ReggeWheeler",
      Resources -> 
     	{
        "Guides/ReggeWheeler"
      }
    }
  }
]

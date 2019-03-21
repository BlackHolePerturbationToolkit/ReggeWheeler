(* ::Package:: *)

BeginPackage["ReggeWheeler`"];

$ReggeWheelerInformation::usage = "$ReggeWheelerInformation is a list of rules that gives information about the version of the ReggeWheeler package you are running.";
$ReggeWheelerInstallationDirectory::usage = "$ReggeWheelerInstallationDirectory gives the top-level directory in which the ReggeWheeler package is installed.";

$ReggeWheelerVersionNumber::usage = "$ReggeWheelerVersionNumber is a real number which gives the current version number for the ReggeWheeler package.";
$ReggeWheelerReleaseNumber::usage = "$ReggeWheelerReleaseNumber is an integer which gives the current release number for the ReggeWheeler package.";
$ReggeWheelerVersion::usage = "$ReggeWheelerVersionNumber is a string that gives the version of the ReggeWheeler package you are running.";

Begin["`Private`"];

$ReggeWheelerInstallationDirectory = FileNameDrop[FindFile["ReggeWheeler`"], -2];

$ReggeWheelerVersionNumber        = 1.0;
$ReggeWheelerReleaseNumber        = 0;

$ReggeWheelerVersion :=
 Module[{path, version, release, buildid, gitrev, gitdir},
  path = $ReggeWheelerInstallationDirectory;
  version = ToString[NumberForm[$ReggeWheelerVersionNumber, {Infinity, 1}]];
  release = ToString[$ReggeWheelerReleaseNumber];

  buildid = Quiet@ReadList[FileNameJoin[{path, "BUILD_ID"}], "String"];
  If[SameQ[buildid, $Failed],
    buildid = "";
  ,
    buildid = " (" <> First[buildid] <> ")";
  ];

  (* First, check for a GIT_REVISION file. If it exists, use its contents as the revision. *)
  gitrev = Quiet@ReadList[FileNameJoin[{path, "GIT_REVISION"}],"String"];

  (* Otherwise, try to determine the git revision directly *)
  If[SameQ[gitrev, $Failed],
    gitdir = FileNameJoin[{path, ".git"}];
    If[FileType[gitdir] === Directory,
      gitrev = Quiet@ReadList["!git --git-dir "<>gitdir<>" rev-parse HEAD", String];
      If[gitrev === {}, gitrev = $Failed];
    ];
  ];

  (* If it worked, ReadList returns a list but we just want the first element (line) *)
  If[Head[gitrev] === List, gitrev = First[gitrev]];

  (* Check we have a git revision and otherwise give up trying *)
  If[Head[gitrev] === String && StringMatchQ[gitrev, RegularExpression["[0-9a-f]{5,40}"]], gitrev = " (" <> gitrev <> ")", gitrev = ""];

  version <> "." <> release <> buildid <> gitrev
]

$ReggeWheelerInformation :=
  {"InstallationDirectory" -> $ReggeWheelerInstallationDirectory,
   "Version" -> $ReggeWheelerVersion,
   "VersionNumber" -> $ReggeWheelerVersionNumber,
   "ReleaseNumber" -> $ReggeWheelerReleaseNumber}

End[];
EndPackage[];

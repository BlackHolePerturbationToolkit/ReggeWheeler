VerificationTest[
  orbit = KerrGeoOrbit[0, 10.`32, 0, 1];
  ReggeWheelerPointParticleMode[2, 2, 2, 0, orbit];
  Names["ReggeWheeler`*`Private`*$" ~~ DigitCharacter ..],
  {},
  TestID -> "Memory leaks"
]

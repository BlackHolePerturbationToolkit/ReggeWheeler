{% include head.html %}

# ReggeWheeler

A Mathematica package for computing solutions to the ReggeWheeler equation. Note this package depends upon the [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/) and [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/) packages to run.

Explicitly the package computes solutions to the Regge-Wheeler equation:

$$ f^2 \frac{d^2\Psi}{dr^2} + f f' \frac{d\Psi}{dr} + (\omega^2 - U_l_(r))\Psi = \mathcal{T} $$

where the potential is given by

$$ U_l(r) = \frac{f}{r^2}\left( l(l+1) -\frac{6M}{r} \right) $$

and
$f = 1-2M/r$  
$\omega$ is the mode frequency  
$\mathcal{T}$ is the source.


Currently the source has been implement for a point particle moving along a circular geodesic orbit. As an example, the flux of gravitational waves in this case for the $l=2,m=2$ mode is easily computed using:

```
r0 = 10`30;
orbit = KerrGeoOrbit[0, r0, 0, 1];

s = 2; l = 2; m = 2; n = 0;
mode = ReggeWheelerPointParticleMode[s, l, m, n, orbit];
```
The `mode` variable is now a `ReggeWheelerMode` object which can be queried for interesting results. For the energy and angular momentum flux use `mode["Fluxes"]` which returns 
```
<|"Energy" -> <|"\[ScriptCapitalI]" -> 0.000026843977395510576377, "\[ScriptCapitalH]" -> 5.6541387345369358451*10^-9|>, "AngularMomentum" -> <|"\[ScriptCapitalI]" -> 0.00084888110027888050954, "\[ScriptCapitalH]" -> 1.78799566077188627708*10^-7|>|>
```

Note the high precision of the input value of $r_0$. This is required by default as the package uses the MST method to calculate the homogeneous solutions.


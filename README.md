# GSW.jl

This Julia implementation of the Thermodynamic Equation of Seawater 2010 (TEOS-10) is based primarily on wrappers of the [GSW-C implementation](https://github.com/TEOS-10/GSW-C).

The original Gibbs-SeaWater (GSW) Oceanographic Toolbox contains the TEOS-10 subroutines for evaluating the thermodynamic properties of pure water (using IAPWS-09) and seawater (using IAPWS-08 for the saline part). Additional information can be obtained  on teos-10.org.

Function names, as well as most argument names, match those used in the [TEOS-10 documentation](http://www.teos-10.org/pubs/gsw/html/gsw_contents.html).

If you use the GSW Oceanographic Toolbox the authors of the original toolbox ask that you include a reference to McDougall and Barker (2011), whose full citation is:
McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5.

## Installing

```
  using Pkg
  Pkg.clone("https://github.com/ax1ine/GSW.jl")
```

## Example
This function calculates in-situ density from Absolute Salinity and Conservative
Temperature, using the computationally-efficient expression for
specific volume in terms of SA, CT and p (Roquet et al., 2014).

```
  using GSW
  rho = gsw_rho(sa, ct, p)
```

where sa - Absolute Salinity (g/kg); ct - Conservative Temperature (ITS-90) (deg C): p - sea pressure (dbar)(i.e. absolute pressure - 10.1325 dbar); rho - in-situ density (kg/m^3)

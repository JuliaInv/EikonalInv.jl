[![Build Status](https://travis-ci.org/JuliaInv/EikonalInv.jl.png?branch=master)](https://travis-ci.org/JuliaInv/EikonalInv.jl) [![Coverage Status](https://coveralls.io/repos/github/JuliaInv/EikonalInv.jl/badge.png?branch=master)](https://coveralls.io/github/JuliaInv/EikonalInv.jl?branch=master)

# EikonalInv.jl
A Julia package for solving the inverse eikonal equation on a regular rectangular mesh.
For forward modelling and sensitivities it uses the fast marching algorithm for the factored eikonal equation.

Based on the following paper (please cite if you are using the package):

Eran Treister and Eldad Haber, A fast marching algorithm for the factored eikonal equation, Journal of Computational Physics, 324, 210-225, 2016.

The package is also used in the following papers in joint inversions:
Lars Ruthotto, Eran Treister and Eldad Haber, jInv--a flexible Julia package for PDE parameter estimation, SIAM Journal on Scientific Computing, 39 (5), S702-S722, 2017. 

Eran Treister and Eldad Haber, Full waveform inversion guided by travel time tomography, SIAM Journal on Scientific Computing, 39 (5), S587-S609, 2017.

# Requirements

This package is intended to use with julia versions 0.4.x.

This package is an add-on for jInv, which needs to be installed.

# Installation

```
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/FactoredEikonalFastMarching.jl","FactoredEikonalFastMarching")
Pkg.clone("https://github.com/JuliaInv/EikonalInv.jl","EikonalInv")

Pkg.test("EikonalInv")
```

# Examples

Under "examples/SEGTravelTimeInversionExample.jl" you can find the 2D experiment that was shown in the paper above. 

# EikonalInv.jl
A Julia package for solving the inverse eikonal equation of a regular rectangular mesh.
For forward modelling and sensitivities it uses the fast marching algorithm for the factored eikonal equation.

Based on the following paper (please cite if you are using the package):

Eran Treister and Eldad Haber, A fast marching algorithm for the factored eikonal equation.

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

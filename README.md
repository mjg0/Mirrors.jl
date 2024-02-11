# Mirrors.jl

This repository contains my thesis and its associated code. It's a Julia package, and can be installed with:

```jldoctest
pkg> add https://github.com/mjg0/Mirrors.jl
```

The thesis, which goes into detail on the research that requires this code, can be built with `make` if you have `pdflatex` and `pygments` installed.

This document will focus more on the code itself since the thesis covers the math and purpose behind `Mirrors.jl`. In short, a `Mirror` is a circular conducting mirror, possibly with some surface roughness:

```julia
using Mirrors, Plots

radius = 5.0 # wavelengths
N = 20 # number of rings; more rings means more precision and more memory use
rms = 0.1 # RMS surface roughness
sigma = 3.0 # standard deviation of surface roughness

M = Mirror(radius, N, rms, sigma)
heatmap(M) # plot the mirror's height
```

The impedance of the mirror (how each point on the mirror interacts with each other point when the mirror is charged) is a matrix of complex doubles, and scales quartically in memory use with the number of rings in the mirror:

```julia
Z = impedance(M)
```

Given the impedance matrix and an electric field on the mirror's surface:

```julia
# Profile of the illuminating beam; default is a uniform plane wave
sigma = 2.0 # wavelengths
function gaussianprofile(r, θ)
    return exp(-r^2/2sigma^2)
end
# Angle of incident beam relative to normal
angle = π/6
# Electric field
E = electricfield(M, angle, beamprofile=gaussianprofile)
heatmap(M, real.(E))
```

...the current on the mirror's surface can be found:

```julia
J = Z \ E
heatmap(M, real.(J))
```

Once the surface current is calculated, reflectance can be found:

```julia
R = Reflectance(M, J)
heatmap(R)
```

A `Reflectance` represents the far field reflectance from a `Mirror` in the mirror's entire upper hemisphere as an [azimuthal equidistant projection](https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection), centered normal to the mirror and extending 90 degrees to the plane of the mirror.
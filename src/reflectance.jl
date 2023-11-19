export reflectancegrid
export reflectancesphere

using SpecialFunctions, HypergeometricFunctions



"""
    reflectanceat(mirror, J, θ, ϕ, λ)

Calculate the far-field reflectance at a certain angle given a mirror and the current thereon

# Arguments
- `mirror::Mirror`: the mirror from which to measure
- `J::Vector{Float64}`: the current on the mirror
- `θ::Real`: the polar angle in radians at which to measure reflectance
- `ϕ::Real`: the azimuthal angle in radians at which to measure reflectance
- `λ::Real`: the wavelength of incident light; default 1
"""
function reflectanceat(mirror::Mirror, J::AbstractVector{ComplexF64}, θ::Real, ϕ::Real, λ::Real=1)
    sum = 0.0+0.0im
    for (cur, (r, t, z)) in zip(J, Iterators.flatten(mirror))
        sum += cur * exp(1im * 2π/λ * (r*cos(ϕ-t)*sin(θ) + z*cos(θ)))
    end
    return sum
end



"""
    expectedreflectanceat(R, α, θ, ϕ, λ)

Calculate expected far-field reflectance at a certain angle for a flat mirror.

If `α` is 0, the simpler Airy formula will be used

# Arguments
- `R::Real`: the radius of the flat mirror
- `θ::Real`: the polar angle in radians at which to measure reflectance
- `ϕ::Real`: the azimuthal angle in radians at which to measure reflectance; default 0
- `α::Real`: the angle of the incident light in radians, measured from normal; default 0
- `λ::Real`: the wavelength of incident light; default 1
"""
function expectedreflectanceat(R::Real, θ::Real=0, ϕ::Real=0, α::Real=0, λ::Real=1)
    # Found this with Mathematica: Integrate[ Exp[I (a r Cos[t] + b r Sin[t])], {r, 0, R}, {t, 0, 2 \[Pi]}, Assumptions -> {{a, b, R} \[Element] Reals, {R, a^2, b^2} > 0}]
    # UPDATE: I missed a multiplication by r to account for the circular integration :facepalm:. When I fixed it, it came out to: \[Pi] R^2 Hypergeometric0F1Regularized[2, 1/4 (a^2 + b^2) R^2]
    if α == 0 # use Airy formula
        c = 2π/λ * R * sin(θ)
        return besselj(1, c) / c
    else
        _₀F₁(a, z) = pFq(typeof(a)[], [a], z)
        return π * R^2 * _₀F₁(2, -(π*R/λ)^2 * (sin(α)^2+sin(θ)^2-2*sin(α)*sin(θ)*cos(ϕ)))
    end
end



defaulta = 0.025
"""
    reflectancegrid(f[, a])
    reflectancegrid(mirror, J[, a])
    reflectancegrid(r[, α[, λ[, a]]])

Calculate the reflectance on a polar azimuthal projection grid that centers at normal and extends to grazing

# Arguments
- `f::Function`: a function of the form `f(θ, ϕ)` (polar angle, azimuthal angle) that returns reflectance at that angle
- `a::Real`: the separation between points on the grid (smaller `a` gives higher resolution); default $defaulta
- `mirror::Mirror`: a `Mirror` for which to calculate the reflectance grid
- `J::Vector{ComplexF64}`: the current on the `Mirror`
- `r::Real`: radius of an ideal mirror for which to calculate theoretical reflectance
- `α::Real`: angle from normal, in radians, of incident beam for ideal mirror; default 0 **WARNING: non-zero α isn't working**
- `λ::Real`: wavelength of incident beam; default 1

# Return Value and Interpretation
`reflectancegrid` returns a 3-tuple containing an x range, a y range, and the reflectance grid itself. The x and y
range are identical and always go from -π/2 to π/2, with `a` determining the step size within each range; both are
returned so that one can use `heatmap(reflectancegrid(f)...)` directly. The grid that is returned is in essence a
circle of radius π/2 (all values outside this circle are 0), with each point on the grid corresponding to the `θ` and
`ϕ` on a hemisphere (which ranges in polar angle from 0 to π/2); in form it resembles a map centered on the north pole
and extending to the equator.

# Examples
Obtain high-resolution expected reflectance for a mirror of radius 15 with incident light angle 30 degrees from normal:
    `x, y, refl = reflectancegrid(15, π/6, 1, 0.005)`

Heatmap the calculated reflectance for a `Mirror` `m` and surface current `J`:
    `heatmap(reflectancegrid(m, J)...)`

Use your own function of (θ, ϕ, λ) to calculate a low-resolution reflectance grid for a wavelength of 2.5:
    `x, y, refl = reflectancegrid((θ, ϕ)->myfunc(θ, ϕ, 2.5), 0.1)`
"""
function reflectancegrid(f::Function, a::Real=defaulta)
    range = -π/2:a:π/2
    return range, range, [sqrt(x^2+y^2)<π/2 ? abs(f(sqrt(x^2+y^2), atan(x, y))) : 0 for x=range, y=range]
end
function reflectancegrid(mirror::Mirror, J::AbstractVector{ComplexF64}, a::Real=defaulta)
    return reflectancegrid((θ, ϕ)->reflectanceat(mirror, J, θ, ϕ), a)
end
function reflectancegrid(r::Real, α::Real=0, λ::Real=1, a::Real=defaulta)
    if α == 0
        return reflectancegrid((θ, ϕ)->airyreflectanceat(r, θ, λ), a)
    else
        return reflectancegrid((θ, ϕ)->expectedreflectanceat(r, α, θ, ϕ, λ), a)
    end
end



"""
    reflectancesphere(mirror, J, λ, a)

Calculate the far-field reflectance of a mirror; returns an array of points of the form (θ, ϕ, reflectance)

# Arguments
- `mirror::Mirror`: the mirror from which to measure
- `J::Vector{Float64}`: the current on the mirror
- `λ::Real`: the wavelength of incident light; default 1
- `a::Real`: the average angular separation between sampled angles; default 0.05

# Examples
Create a scatter plot of reflectance for a mirror `m` and surface current `J` on a hemisphere:
```julia
r = reflectancesphere(m, J)
x = [sin(θ)*cos(ϕ) for (θ, ϕ, refl) in r]
y = [sin(θ)*sin(ϕ) for (θ, ϕ, refl) in r]
z = [cos(θ) for (θ, ϕ, refl) in r]
refl = [refl for (θ, ϕ, refl) in r]
col = ["rgba(\$(elem), 40, 100, 1)" for elem in Int.(round.(128 .+ 127*(abs.(refl)./max(abs.(refl)...))))]
scatter3d(x, y, z, color=col)
```

See https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf for how the points are calculated
"""
function reflectancesphere(mirror::Mirror, J::AbstractVector{ComplexF64}, λ::Real=1, a::Real=0.01)
    θs = a/2:π/round(π/a):π
    ϕs = (θ)->sin(θ)*a/2:2π/round(2π*sin(θ)/a):2π
    return [(θ, ϕ, reflectanceat(mirror, J, θ, ϕ, λ)) for θ in θs for ϕ in ϕs(θ)]
end
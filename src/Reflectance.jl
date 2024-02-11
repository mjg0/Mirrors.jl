export Reflectance

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
    sum = 0.0 + 0.0im
    for (cur, (r, t, z)) in zip(J, Iterators.flatten(mirror))
        sum += cur * exp(1im * 2π / λ * (r * cos(ϕ - t) * sin(θ) + z * cos(θ)))
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
    if α == 0 # use Airy formula
        c = 2π / λ * R * sin(θ)
        return besselj(1, c) / c
    else
        _₀F₁(a, z) = pFq(typeof(a)[], [a], z)
        return π * R^2 * _₀F₁(2, -(π * R / λ)^2 * (sin(α)^2 + sin(θ)^2 - 2 * sin(α) * sin(θ) * cos(ϕ)))
    end
end



defaultn = 200
"""
    Reflectance

A struct containing a polar azimuthal reflectance grid that extends to grazing.

`Reflectance` is an `AbstractMatrix` containing the reflectance grid, which extends from
`-π/2` to `π/2`; it resembles a map of the earth that extends to the equator, centered at
the north pole. `heatmap` and `write` are overloaded for `Reflectance`.

# Examples

Obtain high-resolution expected reflectance for a mirror of radius 15 with incident light angle 30 degrees from normal:
    `refl = Reflectance(15, π/6, 1, 1000)`

Heatmap the calculated reflectance for a `Mirror` `m` and surface current `J`:
    `heatmap(Reflectance(m, J))`

Use your own function of (θ, ϕ, λ) to calculate a low-resolution reflectance grid for a wavelength of 2.5:
    `refl = Reflectance((θ, ϕ)->myfunc(θ, ϕ, 2.5), 0.1)`

---

    Reflectance(f[, a])
    Reflectance(mirror, J[, a])
    Reflectance(r[, α[, λ[, a]]])
    Reflectance(io)
    Reflectance(filename)

# Arguments
- `f::Function`: a function of the form `f(polar, azimuthal)` returing reflectance there
- `n::Integer`: the number of points along each axis of the grid; default $defaultn
- `mirror::Mirror`: a `Mirror` for which to calculate the reflectance grid
- `J::Vector{ComplexF64}`: the current on the `Mirror`
- `r::Real`: radius of an ideal mirror for which to calculate theoretical reflectance
- `α::Real`: angle from normal, in radians, of incident beam for ideal mirror; default 0
- `λ::Real`: wavelength of incident beam; default 1
- `io::IO`: an `IO` to which to write the `Reflectance`
- `filename::AbstractString`: a file to which to write the `Reflectance`
"""
struct Reflectance <: AbstractMatrix{Float64}
    r::Matrix{Float64}

    function Reflectance(f::Function, n::Integer=defaultn)
        grid = range(-π / 2, π / 2, length=n)
        return new([sqrt(x^2 + y^2) < π / 2 ? abs(f(sqrt(x^2 + y^2), atan(x, y))) : 0
                    for x in grid, y in grid])
    end

    function Reflectance(mirror::Mirror, J::AbstractVector{ComplexF64}, n::Integer=defaultn)
        return Reflectance((θ, ϕ) -> reflectanceat(mirror, J, θ, ϕ), n)
    end

    function Reflectance(r::Real, α::Real=0, λ::Real=1, n::Integer=defaultn)
        return Reflectance((θ, ϕ) -> expectedreflectanceat(r, θ, ϕ, α, λ), n)
    end

    function Reflectance(io::IO)
        n = read(io, UInt64)
        r = Matrix{Float64}(undef, n, n)
        read!(io, r)
        return new(r)
    end

    Reflectance(filename::AbstractString) = Reflectance(open(filename))
end



# Implement AbstractVector interface
Base.size(refl::Reflectance) = size(refl.r)

Base.getindex(refl::Reflectance, args...) = getindex(refl.r, args...)

Base.setindex!(refl::Reflectance, args...) = setindex!(refl.r, args...)



# Write a Reflectance to an IO
function Base.write(io::IO, refl::Reflectance)
    write(io, UInt64(first(size(refl))))
    write(io, refl.r)
end
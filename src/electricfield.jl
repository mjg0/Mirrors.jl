export electricfield

using LinearAlgebra



"""
    electricfield(mirror::Mirror, α, λ=1; beamprofile=(r, θ)->1.0)

Return the electric field on `mirror` given an incident angle `α` (in radians, measured from
normal) and wavelength `λ`.

`beamprofile` is a function taking two arguments, a distance and an angle; this defines a
point relative to the center of the beam. The output of the function is multiplied by the
natural reflectance of the point on the surface that the given part of the beam would
strike. If `beamprofile` is not specified, it defaults to a function producing a plane wave.

Ordering of the electric field points is the same as that of `Iterators.flatten(mirror)`.

See also [`impedance`](@ref), [`surfacecurrent`](@ref).
"""
function electricfield(mirror::Mirror, α::Real, λ::Real=1.0;
                       beamprofile::Function=(r, θ)->1.0)
    kx::T =  2π/λ*sin(α)
    kz::T = -2π/λ*cos(α)
    return [begin
        X = r*cos(θ)*cos(α)-z*sin(α)
        Y = r*sin(θ)
        beamprofile(norm((X, Y)), atan(Y/X)) * exp(im*(kx*r*cos(θ)+kz*z))
    end for (r, θ, z) in Iterators.flatten(mirror)]
end

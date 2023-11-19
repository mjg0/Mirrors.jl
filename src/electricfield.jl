export electricfield



"""
    electricfield(mirror, α, λ)

Return the electric field on a mirror given an incident angle and wavelength

# Arguments
- `mirror::Mirror`: the mirror on which to calculate the electric field
- `α::Real`: the incident angle, measured from normal, in radians
- `λ::Real`: the wavelength of the incident light
"""
function electricfield(mirror::Mirror, α::Real, λ::Real=1)
    kx =  2π/λ * sin(α)
    kz = -2π/λ * cos(α)
    return [exp(im * (kx * r * cos(θ) + kz * z)) for (r, θ, z) in Iterators.flatten(mirror)]
end

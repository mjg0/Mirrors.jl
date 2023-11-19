using StaticArrays



"""
Generator that returns all `(m, n)` (annulus and patch) indices consituting a mirror for a given number of `rings`
"""
function mirrorindices(rings)
    return ((m, n) for n=0:rings-1 for m=0:6*n+2)
end



"""
Given annular width `a` and ring index `n`, return the radii of the outer and inner points of a patch
"""
function rplusminus(a::Real, n::Integer)
    r = a * (n + 0.5)
    dr = a / 2sqrt(3)
    return r + dr, r - dr
end



"""
Given patch index `m` and ring index `n`, return the two angles of the points of a patch
"""
function θplusminus(m::Integer, n::Integer)
    θ = 2π * (m + 0.5) / (6 * n + 3)
    dθ = π / sqrt(3) / (6 * n + 3)
    return θ + dθ, θ - dθ
end



@doc raw"""
`Patch`

A struct consisting of the radius, angle, height, and surface Jacobian of the 4 points of a patch

Patches are pieces of a circle, each having equal area. One patch is defined by an inner and outer radius and two
angles. A patch is represented by 4 points (2 radii and 2 angles, somewhat in from the edges of each patch), each with a
radius, angle, height, and surface Jacobian. The radius and angle of each of these points is determined by the annular
width of each ring, which ring (`n`) the patch is in, and which patch on the ring (`m`) the points are in.

A `Patch` constitutes 4 `SVector`s, `r` (radius), `θ` (angle), `z` (height), and `s` (surface Jacobian), each with 4
elements corresponding to each of the 4 points.

Given middle radius and angle (`r` and `θ`) and annular and angular widths (`a` and `w`), the 4 points are (in order):

```math
(r-\frac{a}{2\sqrt 3}, \theta-\frac{w}{2\sqrt 3})
(r+\frac{a}{2\sqrt 3}, \theta-\frac{w}{2\sqrt 3})
(r+\frac{a}{2\sqrt 3}, \theta+\frac{w}{2\sqrt 3})
(r-\frac{a}{2\sqrt 3}, \theta+\frac{w}{2\sqrt 3})
```

This means that points on the patch are indexed thus:

```
   ----------
   | 2    3 |
   |        |
^  | 1    4 |
r  ----------
   θ >
```
"""
struct Patch <: AbstractVector{NTuple{3,Float64}}
    r::SVector{4,Float64}
    θ::SVector{4,Float64}
    z::SVector{4,Float64}
    m::UInt32
    n::UInt32
    """
    `Patch(a, m, n, z)`

    Construct a `Patch` corresponding to the patch determined by annular width `a`, patch `m`, and ring `n`

    Arguments:
    - `a`: annular width of each ring
    - `m`: which patch on the ring this will be
    - `n`: which ring this patch will be in
    - `z`: a function, taking a radius and angle, that determines the height at that point
    """
    function Patch(a::Real, m::Integer, n::Integer, z::Function=(r,θ)->0)
        rplus, rminus = rplusminus(a, n)
        θplus, θminus = θplusminus(m, n)
        r  = @SVector [rminus,  rplus,   rplus,  rminus]
        θ  = @SVector [θminus, θminus, θplus, θplus]
        z_ = @SVector [z(r[1], θ[1]), z(r[2], θ[2]), z(r[3], θ[3]), z(r[4], θ[4])]
        return new(r, θ, z_, m, n)
    end
end



# Implement abstract vector interface for patch
Base.size(::Patch) = Tuple(4)

Base.getindex(patch::Patch, i) = patch.r[i], patch.θ[i], patch.z[i]

function Base.setindex!(patch::Patch, v::NTuple{3,<:Real}, i)
    patch.r[i] = v[1]
    patch.θ[i] = v[2]
    patch.z[i] = v[3]
end

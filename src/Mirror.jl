export Mirror

using Statistics, Interpolations, FFTW, Serialization



@doc raw"""
    noisy2dspline(len, n, rms, σ)

Return a square CubicSplineInterpolation with side `len`, `n` points, RMS height `rms`, roughness determined by `σ`

Since roughness on real surfaces tends to fall off as a Gaussian at higher frequencies, an array of random numbers is
generated, put through a Gaussian low-pass filter, and returned as a CubicSplineInterpolation. `σ` is the standard
deviation--higher `σ` means higher-frequency roughness.

- `len`: the side length of the square surface
- `n`: dimension of the grid of points, which will have $n^2$ total points
- `rms`: RMS height of the surface
- `σ`: standard deviation of roughness in frequency space
"""
function noisy2dspline(len::Real, n::Integer, rms::Real, σ::Real)
    # grid of points
    x = range(-len/2, len/2, length=n)
    y = range(-len/2, len/2, length=n)
    # frequency-space grid
    dk = 2π / n
    k1d = 0.5*dk-π:dk:π # circshift(0.5*dk-π:dk:π, ceil(n/2))
    k = k1d * reshape(k1d, (1, n))
    # Fourier z
    zf = fft(randn(n, n), (1, 2)) .* exp.(-k.^2/2/(σ*dk)^2)
    # real component of inverse transform of `zf` is the collection of points
    z = real(ifft(zf))
    # return the resultant spline interpolation, correcting RMS and zeroing the average surface height
    offset = mean(z)
    zrms = sqrt(sum((z.-offset).^2)) / n
    return cubic_spline_interpolation((x, y), (z.-offset).*(rms/zrms))
end



"""
    Mirror

A struct containing the array of `Patch`s that make up a mirror.

A `Mirror` constitutes the patch annular width `a`, number of rings `rings`, and `Patch` array `patches`.

A `Mirror` represents a circular conducting surface, possibly with some height, that is split into patches of equal
area. The first ring of the mirror (the middle) is simply the inner circle of radius `a` and has 3 patches, resembling a
pie chart with 3 equal areas. The next is a ring of inner radius `a` and outer radius `2a`, which is split into 9
patches of equal angular width. The next ring has 15 patches, the next 21, and so on.

Since storing the array of `Patch`s is the main concern of `Mirror`, the `iterator` and `getindex` operators are
overloaded for convenience:

```julia-repl
julia> m = Mirror(1, 1)
julia> typeof(m[1])
Mirrors.Patch
julia> for p in m
           println(p.z)
       end
[0.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, 0.0]
```
"""
struct Mirror <: AbstractVector{Patch}
    a::Float64
    rings::Int64
    patches::Vector{Patch}
    z::Function
    s::Function
    """
        Mirror(r, rings, z, s)

    Construct a `Mirror` of radius `r` with `rings` total rings, with `z` height and `s` surface Jacobian

    Arguments:
    - `rings`: number of rings on the mirror
    - `r`: radius of the mirror
    - `z`: function, taking a radius and an angle, that determines the height at that point
    - `s`: function, taking a radius and angle, that determines the surface Jacobian at that point
    """
    function Mirror(r::Real, rings::Integer, z::Function, s::Function)
        a = r / rings
        patches = [Patch(a, m, n, z) for (m, n) in mirrorindices(rings)]
        return new(a, rings, patches, z, s)
    end # function Mirror
    """
        Mirror(r, rings, rms, σ)

    Construct a `Mirror` of radius `r` with `rings` total rings, with roughness defined by `rms` and `σ`

    `rms` represents the RMS height of the rough mirror. `σ` is a measure of the frequency of the mirror; higher
    `σ` means higher-frequency roughness. The mirror roughness is essentially determined by generating many random
    points then putting then through a low-pass filter--see `noisy2dspline`

    Arguments:
    - `rings`: number of rings on the mirror
    - `r`: radius of the mirror
    - `rms`: RMS height of the mirror; default 0
    - `σ`: standard deviation of roughness in frequency space; default 0
    """
    function Mirror(r::Real, rings::Integer, rms::Real=0, σ::Real=0)
        z = (r, θ) -> 0
        s = (r, θ) -> 1
        if rms != 0
            spline = noisy2dspline(2*r, 8*rings, rms, σ)
            z = (r::Real, θ::Real) -> spline(r*cos(θ), r*sin(θ))
            s = (r::Real, θ::Real) -> sqrt(1 + sum(Interpolations.gradient(spline, r*cos(θ), r*sin(θ)).^2))
        end
        return Mirror(r, rings, z, s)
    end # function Mirror
    """
        Mirror(io)
        Mirror(filename)

    Construct a `Mirror` from a serialized file or stream
    """
    Mirror(io::IO) = deserialize(io)
    Mirror(filename::AbstractString) = deserialize(filename)
end # struct Mirror



# Implement AbstractVector interface
Base.size(mirror::Mirror) = size(mirror.patches)

Base.getindex(mirror::Mirror, args...) = getindex(mirror.patches, args...)

Base.setindex!(mirror::Mirror, v::Patch, args...) = setindex!(mirror.patches, v, args...)



# Write a mirror to a stream or file
Base.write(io::IO, mirror::Mirror) = serialize(io, mirror)

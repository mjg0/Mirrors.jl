export Mirror

using Statistics, Interpolations, FFTW, Serialization, ImageFiltering



# Magic bytes for special mirrors
roughmirrorcode::UInt16 = 0x01



"""
    roughmirrorz_h(h, r)

Return `z` and `s` functions given height array `h` and mirror radius `r`.
"""
function roughmirror_z_s(h::AbstractArray, r::Real)
    span = range(-r, r, length=first(size(h)))
    itp = scale(interpolate(h, BSpline(Cubic(Free(OnCell())))), span, span)
    z = (r, θ) -> itp(r*cos(θ), r*sin(θ))
    s = (r, θ) -> sqrt(1 + sum(Interpolations.gradient(itp, r*cos(θ), r*sin(θ)).^2))
    return z, s
end



"""
    Mirror

A struct containing the array of `Patch`s that make up a mirror.

A `Mirror` constitutes the patch annular width `a`, number of rings `rings`, and `Patch` array `patches`.

A `Mirror` represents a circular conducting surface, possibly with some height, that is split into patches of equal
area. The first ring of the mirror (the middle) is simply the inner circle of radius `a` and has 3 patches, resembling a
pie chart with 3 equal areas. The next is a ring of inner radius `a` and outer radius `2a`, which is split into 9
patches of equal angular width. The next ring has 15 patches, the next 21, and so on.

Mirrors can be serialized to and from files and `IO`s using the constructor and `write`.

Since storing the array of `Patch`s is the main concern of `Mirror`, it is an `AbstractArray{Patch}`:

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
            # Interpolate with 8*rings points along each axis
            n = 8*rings
            # Make a supergrid and subgrid to prevent sharp edges
            h₀ = imfilter(rand(2n, 2n).-0.5, Kernel.gaussian(σ*n/2r))
            h = h₀[1+n÷2:2n-n÷2,1+n÷2:2n-n÷2]
            # Transform Z appropriately
            h .-= mean(h)
            h .*= rms/sqrt(sum(x->x^2, h))*n
            # z and s given h and r
            z, s = roughmirror_z_s(h, r)
        end
        return Mirror(r, rings, z, s)
    end # function Mirror
    """
        Mirror(io)
        Mirror(filename)

    Construct a `Mirror` from a file or stream
    """
    function Mirror(io::IO)
        a = read(io, Float64)
        rings = read(io, Int64)
        patches = Vector{Patch}(undef, length(collect(mirrorindices(rings))))
        read!(io, patches)
        if eof(io)
            return new(a, rings, patches, (args...)->0, (args...)->1)
        end
        magicbytes = read(io, UInt16)
        if magicbytes == roughmirrorcode
            n = read(io, Int64)
            height = Matrix{Float64}(undef, n, n)
            read!(io, height)
            z, s = roughmirror_z_s(height, a*rings)
            return new(a, rings, patches, z, s)
        else
            throw(ErrorException("Mirror appears corrupt"))
        end
    end
    Mirror(filename::AbstractString) = Mirror(open(filename))
end # struct Mirror



# Implement AbstractVector interface
Base.size(mirror::Mirror) = size(mirror.patches)

Base.getindex(mirror::Mirror, args...) = getindex(mirror.patches, args...)

Base.setindex!(mirror::Mirror, v::Patch, args...) = setindex!(mirror.patches, v, args...)



# Write a mirror to a stream or file
function Base.write(io::IO, mirror::Mirror)
    byteswritten = (write(io, mirror.a)
                   +write(io, mirror.rings)
                   +write(io, mirror.patches))
    if hasproperty(mirror.z, :itp)
        byteswritten += (write(io, roughmirrorcode)
                        +write(io, first(size(mirror.z.itp)))
                        +write(io, [z for z in mirror.z.itp]))
    end
    return byteswritten
end

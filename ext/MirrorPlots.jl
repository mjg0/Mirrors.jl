module MirrorPlots



using Plots



"""
Return the point (and its index) on `mirror` closest to `(r, th)`

The point is of the form `(r, θ, z)`, and the index corresponds to which point on the mirror it is.

For example, `closestpoint(m, 1, 3)` might return `((1.1, 2.5, -1.4), 6)` (the 2nd point on the 2nd patch)
"""
function closestpoint(mirror::Mirror, r::Real, θ::Real)
    n = Int(round(min(r / mirror.a - 0.5, mirror.rings-1)))
    m = Int(round(mod(θ, 2π) / 2π * (6n + 3) - 0.5))
    i = findall(idx->idx==(m,n), [idx for idx in mirrorindices(mirror.rings)])[1]
    patch = mirror[i]
    pointdists = [sqrt(rp^2+r^2-rp*r*2cos(θp-θ)) for (rp,θp,zp) in patch]
    j = argmin(pointdists)
    return patch[j], 4i-4+j
end



# Overload Plots.heatmap for Mirror
function Plots.heatmap(mirror::Mirror, height::Union{AbstractVector{<:Real},Nothing}=nothing#;
                       #resolution=200, color=:bluesreds, clims=nothing
                       , kw...)
    height = height === nothing ? [z for (r,θ,z) in Iterators.flatten(mirror)] : height
    r = mirror.rings * mirror.a
    xs = ys = range(-r, r, length=resolution)
    z = [x^2+y^2<r^2 ? height[closestpoint(mirror, sqrt(x^2+y^2), atan(y, x))[2]] : 0 for x in xs, y in ys]
    minmax = max(abs.(z)...)
    clims = clims === nothing ? (-minmax-0.001, minmax+0.001) : clims
    heatmap(xs, ys, z; color=color, clims=clims, kw...)
end



end # module

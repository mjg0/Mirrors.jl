export impedance

using LinearAlgebra, HCubature



"""
    greens(r1, θ1, z1, r2, θ2, z2, z2[, λ=1.0])

The 3-dimensional Greens function using cylindrical coordinates given wavelength `λ`.
"""
function greens(r1::Real, θ1::Real, z1::Real, r2::Real, θ2::Real, z2::Real, λ::Real=1.0)
    x1, y1 = r1.*(cos(θ1), sin(θ1))
    x2, y2 = r2.*(cos(θ2), sin(θ2))
    ρ = norm((x1, y1, z1) .- (x2, y2, z2))
    return exp(2π/λ * 1im * ρ) / (4π * ρ)
end



"""
    TransformParameters{A, S, T}

Store the B matrix parameters `a`, `s`, and `t` in a type.
"""
struct TransformParameters{A, S, T} end



"""
    transformfunction(f, zfunc, sfunc, u, v, m, n, tp::TransformParameters)

Return `u*r * f(r, θ, z(r, θ)) * s(r, θ)*J`, transforming `u` and `v` to `r` and `θ` as
defined by `m`, `n`, and `tp`.
"""
@generated function uvtransform(f, zfunc, sfunc, u, v, m, n,
                                ::TransformParameters{A, S, T}) where {A, S, T}
    # B matrix
    poffset = S == 1 || S == 2 ? 1 : -1
    qoffset = S == 1 || S == 4 ? 1 : -1
    p2 = poffset + (T == 3 || T == 4 ? √3 : -√3)
    p3 = poffset + (T == 2 || T == 3 ? √3 : -√3)
    q2 = qoffset + (T == 2 || T == 3 ? √3 : -√3)
    q3 = qoffset + (T == 1 || T == 2 ? √3 : -√3)
    B11, B12, B21, B22 = (p2, (p3-p2),
                          q2, (q3-q2)) .* A ./ 2√3
    # Jacobian (needs to be divided by 6n+3)
    J0 = 2π/A*abs(det([B11 B12;B21 B22]))
    # r and θ offsets
    roffset = S == 2 || S == 3 ? A/2√3 : -A/2√3
    θoffset = S == 3 || S == 4 ? A/2√3 : -A/2√3
    return quote
        B = @SMatrix [$B11 $B12
                      $B21 $B22]
        J = $J0/(6n+3)
        r = B[2,1]*u + B[2,2]*u*v + $roffset + A*(n+0.5)
        θ = 2π * (A*(m+0.5) + B[1,1]*u + B[1,2]*u*v + $θoffset) / (A*(6n+3))
        return u * r * f(r, θ, zfunc(r, θ)) * sfunc(r, θ) * J
    end
end



"""
    integratepatch(f::Function, mirror, patch, s)

Return the integral of `f(r, θ, z)` over the `s`th point on `patch`.
"""
function integratepatch(f, mirror, patch, s)
    total = 0.0
    zf::Function = mirror.z
    sf::Function = mirror.s
    for t=1:4  # total of all 4 triangles
        tp = TransformParameters{mirror.a, s, t}()
        total += hcubature(uv->uvtransform(f, zf, sf, uv[1], uv[2], patch.m, patch.n, tp),
                           (0.0, 0.0), (1.0, 1.0))[begin]
    end
    return total
end



"""
    singularweightelements(f, mirror, patch, s)

Return the weights given function `f(r1, θ1, z1, r2, θ2, z2)` for the `s`th row of the
singular block of the impedance matrix corresponding to `patch` as a tuple.
"""
function singularweightelements(f, mirror, patch, s)
    a = mirror.a
    p = patch[s]
    # K integrals: integrating f(), x*f(), y*f(), and x*y*f()
    K   = integratepatch((r, θ, z)->f(p..., r, θ, z),                   mirror, patch, s)
    Kx  = integratepatch((r, θ, z)->f(p..., r, θ, z)*r*cos(θ),          mirror, patch, s)
    Ky  = integratepatch((r, θ, z)->f(p..., r, θ, z)*r*sin(θ),          mirror, patch, s)
    Kxy = integratepatch((r, θ, z)->f(p..., r, θ, z)*r^2*cos(θ)*sin(θ), mirror, patch, s)
    # return the weight elements
    return K/4, sqrt(3)*Kx/2a, sqrt(3)*Ky/2a, 3Kxy/a^2
end



"""
    singularblockfill!(f, Zblock, mirror, patch)

Fill a 4x4 block of the main diagonal of `Z` given function `f(r1, θ1, z1, r2, θ2, z2)`.
"""
function singularblockfill!(f, Zblock::AbstractMatrix{ComplexF64}, mirror, patch)
    # loop over points in this patch
    @inbounds for (s, p) in enumerate(patch)
        # elements that will make up weights
        e1, e2, e3, e4 = singularweightelements(f, mirror, patch, s)
        # put weights into this strip of Z
        Zblock[s,1] = e1 - e2 - e3 + e4
        Zblock[s,2] = e1 - e2 + e3 - e4
        Zblock[s,3] = e1 + e2 + e3 + e4
        Zblock[s,4] = e1 + e2 - e3 - e4
    end
end



"""
    nonsingularblockfill!(f, Zblock, mirror, patch1, patch2)

Fill a 4x4 block off the main diagonal of `Z` given function `f(r1, θ1, z1, r2, θ2, z2)`.
"""
function nonsingularblockfill!(f, Zblock::AbstractMatrix{ComplexF64}, mirror, patch1, patch2)
    J = π * mirror.a / (12 * patch2.n + 6) # Jacobian
    # 2d loop over points in each patch
    @inbounds for (i, p1) in enumerate(patch1), (j, p2) in enumerate(patch2)
        # r2 * f(p1, p2) * s(p2) * J
        Zblock[i,j] = p2[1] * f(p1..., p2...) * mirror.s(p2[1], p2[2]) * J
    end
end



"""
    pidx(i)

Return the patch index (∈[1,npatches] and the point index (∈[1,4]) within that patch of the `i`th point on a `Mirror`

# Examples
`pidx(7)` returns `(2, 3)` since the 7th point on the `Mirror` is the third point of the second patch.

To get the point corresponding to `i` you can use `mirror[pidx(i)...]`.
"""
function pidx(i::Integer)
    return (i+3)>>2, 1+(i-1)%4
end



"""
    impedance(mirror[, λ=1, singular=true])
    impedance(mirror, f, [singular=true])

Return the impedance matrix corresponding to `mirror` and `f(r1, θ1, z1, r2, θ2, z2)`.

If `λ` is supplied, the function will be `greens(..., λ)`; by default it's `greens(..., 1)`.
"""
function impedance(mirror::Mirror, f::Function=greens, singular::Bool=true)
    # allocate Z
    n = 4 * length(mirror)
    Z = Matrix{ComplexF64}(undef, n, n)
    # fill Z; have to use `collect` because @threads isn't mature yet :/
    @inbounds @fastmath Threads.@threads for (i, patch1) in collect(enumerate(mirror))
        for (j, patch2) in enumerate(mirror)
            Zblock = view(Z, 4i-3:4i, 4j-3:4j)
            if singular && i==j # main diagonal
                singularblockfill!(f, Zblock, mirror, patch1)
            else                # off-diagonal
                nonsingularblockfill!(f, Zblock, mirror, patch1, patch2)
            end
        end
    end
    return Z
end

function impedance(mirror::Mirror, λ::Real, singular::Bool=true)
    return impedance(mirror, (args...)->greens(args..., λ), singular)
end

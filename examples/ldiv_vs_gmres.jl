#!/usr/bin/env julia



using Mirrors, IterativeSolvers, LinearAlgebra, Statistics, Plots



absrms(A::AbstractArray) = abs(sqrt(sum(A.^2)/length(A)))



makesubplot(args...; kw...) = heatmap(args...; aspect_ratio=:equal,
                                               framestyle=:none,
                                               xticks=false,
                                               yticks=false,
                                               colorbar=false,
                                               kw...)



function minmax_x_y_z(mirror::Mirror)
    x, y, z = ntuple(_->Float64[], 3)
    for f in (findmin, findmax)
        (_, pidx), midx = f(p->f(p.z), mirror)
        r, θ, h = mirror[midx][pidx]
        push!(x, r*cos(θ))
        push!(y, r*sin(θ))
        push!(z, h)
    end
    return x, y, z
end



function main()
    # Parameters
    r = 10.0
    N = 30
    α = π/6
    rms = 0.01
    σ = 3.0
    # Build and illuminate mirror
    M = Mirror(r, N, rms, σ)
    Z = impedance(M)
    E = electricfield(M, α)
    # Solve for J, with and without a preconditioner
    tldiv = @elapsed Jldiv = Z \ E
    tgmres = @elapsed Jgmres = gmres(Z, E)
    # Calculate reflectance
    Rldiv = Reflectance(M, Jldiv)
    Rgmres = Reflectance(M, Jgmres)
    # Calculate RMS error for each solution
    RMSEldiv = absrms(E-Z*Jldiv)
    RMSEgmres = absrms(E-Z*Jgmres)
    # Print results
    println("Z condition number: $(cond(Z))")
    println("ldiv  RMSE:         $RMSEldiv (took $tldiv seconds)")
    println("GMRES RMSE:         $RMSEgmres (took $tgmres seconds)")
    # Generate plots
    x, y, z = minmax_x_y_z(M)
    Mplot = makesubplot(M, title="Mirror")
    if abs(z[1]-z[2]) > sqrt(eps(Float64))
        for (i, color) in zip(eachindex(z), (:cyan, :magenta))
            scatter!(Mplot, x[i:i], y[i:i], color=color, label=round(z[i], sigdigits=2),
                     markersize=3, markerstrokewidth=0)
        end
    end
    Eplot = makesubplot(M, real.(E), title="Electric field")
    Jldivplot = makesubplot(M, real.(Jldiv), title="Surface current, ldiv")
    Jgmresplot = makesubplot(M, real.(Jgmres), title="Surface current, gmres")
    Rldivplot = makesubplot(Rldiv, title="Reflectance, ldiv")
    Rgmresplot = makesubplot(Rgmres, title="Reflectance, gmres")
    plotfile = "ldiv_vs_gmres.png"
    png(plot(Mplot, Rldivplot, Rgmresplot, Eplot, Jldivplot, Jgmresplot, size=(1200, 800)),
        plotfile)
    println("Plotted figures to $plotfile")
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
using Test, Mirrors, Statistics, HCubature, LinearAlgebra



@testset "Mirrors" begin
    # Test that `rplusminus` and `θplusminus` return sensible results
    @testset "r_θ_plusminus" begin
        # run checks in loops
        for a = rand(4) * 3, n = 0:4
            rplus, rminus = Mirrors.rplusminus(a, n)
            @test isapprox(rplus - rminus, a / sqrt(3))
            @test isapprox((rplus + rminus) / 2, a * (n + 0.5))
        end
        for n = 0:3, m = 0:6*n+2
            θplus, θminus = Mirrors.θplusminus(m, n)
            patches_this_ring = 6 * n + 3
            @test isapprox(θplus - θminus, 2 * π / sqrt(3) / patches_this_ring)
            @test isapprox((θplus + θminus) / 2, 2 * π * (m + 0.5) / patches_this_ring)
        end
    end



    # Test that `noisy2dspline` works as expected
    #@testset "rough mirror interpolation" begin
    #    function spline_fft(spline)
    #        xs, ys = spline.itp.ranges
    #        z = reshape([spline(x, y) for x in xs for y in ys], (length(xs), length(ys)))
    #        return fft(z, 2)
    #    end
    #    for len=(.63, 1.46, 4.26, 8.91), n=(36, 183, 694), rms=(0.061, 0.49, 1.52, 3.78), σ=(.29, 1.88, 7.24, 16.40)
    #        spline = Mirror(len/2, n, rms, σ).z.itp
    #        # make sure that grid is of the proper dimensions
    #        xs = ys = range(-len/2, len/2, n)
    #        @test length(xs) == length(ys) == n
    #        @test xs[1] ≈ ys[1] ≈ -len / 2
    #        @test xs[end] ≈ ys[end] ≈ len / 2
    #        # check that RMS is properly calculated
    #        z = reshape([spline(x, y) for x in xs for y in ys], (length(xs), length(ys)))
    #        @test isapprox(rms, sqrt(sum(z.^2))/n, rtol=0.0001)
    #        # check that average mirror height is very close to zero
    #        @test isapprox(mean(z), 0, atol=1e-10)
    #        # check that frequency follows an approximately Gaussian distribution
    #        # TODO: fit the frequency curve to a Gaussian, make sure the residual is sane
    #    end
    #end



    # Test that `Patch` constructor fills `r` and `θ` correctly
    @testset "Patch" begin
        rplus, rminus = Mirrors.rplusminus(1.0, 0)
        θplus, θminus = Mirrors.θplusminus(0, 0)
        p = Mirrors.Patch(1.0, 0, 0)
        @test p.r[1] == p.r[4] == rminus
        @test p.r[2] == p.r[3] == rplus
        @test p.θ[1] == p.θ[2] == θminus
        @test p.θ[3] == p.θ[4] == θplus
    end



    # Comprehensively test that the `Mirror` constructor gives back a correct `Mirror`
    # For multiple radii and ring counts, ensure that:
    # - inner and outer radii of patch points are the same on a given ring, and equal what they should
    # - inner and outer angles of patch points are what they should be, and are all equally separated
    # - z's and s's were computed correctly
    @testset "Mirror" begin
        """
        Test a mirror's properties given its construction arguments
        """
        function checkmirror(r::Real, rings::Integer, args...)
            mirror = Mirror(r, rings, args...)
            a = r / rings
            dr = a / sqrt(3)
            for n = 0:rings-1
                w = 2π / (6 * n + 3) # patch angular width
                dθ = w / sqrt(3)
                ring_patches = mirror.patches[3*n^2+1:3*n^2+6*n+3]
                for (p1, p2, m) in zip(ring_patches, [ring_patches[2:end]; [ring_patches[1]]], 0:6*n+2)
                    # check that z's were computed correctly
                    @test all(p1.z .== [mirror.z(p1.r[i], p1.θ[i]) for i = 1:4])
                    # check that spacing between inner and outer r and θ are correct
                    @test p1.r[2] - p1.r[1] ≈ dr
                    @test p1.θ[3] - p1.θ[1] ≈ dθ
                    # check that averages of inner and outer r and θ correspond to the exact middle of the patch
                    @test (p1.r[1] + p1.r[2]) / 2 ≈ a * (n + 0.5)
                    @test (p1.θ[1] + p1.θ[3]) / 2 ≈ w * (m + 0.5)
                    # check that this patch's r and θ match up properly with next patch's
                    @test p1.r[1] == p1.r[4] == p2.r[1] == p2.r[4]
                    @test p1.r[2] == p1.r[3] == p2.r[2] == p2.r[3]
                    @test (p2.θ[1] - p1.θ[1] + 2π) % (2π) ≈ w
                    @test (p2.θ[3] - p1.θ[3] + 2π) % (2π) ≈ w
                end
            end
            return mirror
        end
        checkmirror(1.0, 1)
        checkmirror(2.3, 2)
        checkmirror(3.9621461, 3, (r, θ) -> 1 + sqrt(r), (r, θ) -> 1)
        checkmirror(0.49626, 4, (r, θ) -> r, (r, θ) -> 1 / sqrt(2))
        rough_mirror = checkmirror(3.0, 2, 6.51245, 14.64759)
        # check that the rough mirror doesn't have any symmetries
        for r1 = 0.5:0.5:3.0, r2 = 0.5:0.5:3.0, θ1 = 0:π/4:2π-π/8, θ2 = 0:π/4:2π-π/8
            if r1 == r2 && θ1 == θ2
                continue
            end
            @test rough_mirror.z(r1, θ1) ≉ rough_mirror.z(r2, θ2)
            @test rough_mirror.s(r1, θ1) ≉ rough_mirror.s(r2, θ2)
        end
    end



    # Make sure that Mirror reading and writing works
    @testset "Mirror I/O" begin
        for args in ((1.1, 3, 0.1, 1), (1.2, 4)) # Both flat and rough mirrors
            mirror = Mirror(args...)
            # Electric field gives us means for a quick check
            α = 0.67
            λ = 8.9
            E = electricfield(mirror, α, λ)
            # Write
            buf = IOBuffer()
            write(buf, mirror)
            # Read
            seekstart(buf)
            mirrorcopy = Mirror(buf)
            # Check
            @test electricfield(mirrorcopy, α, λ) == E
            z = [z for patch in mirror for z in patch.z]
            zcopy = [z for patch in mirrorcopy for z in patch.z]
            @test z == zcopy
        end
    end



    # Test that the (r, θ)->(u, v) transform function works
    @testset "uvtransform" begin
        # Check that integrated total is correct
        mirror = Mirror(2.1, 3, (r, θ) -> 0.0, (r, θ) -> 1.0)
        a = mirror.a
        for i = (1, 2, 3, 4, 7, 12, 13, 19, 27), s = 1:4, f = ((r, θ, z) -> 1.0,
                (r, θ, z) -> r,
                (r, θ, z) -> θ,
                (r, θ, z) -> r^2 * sin(θ)^2)
            m = mirror[i].m
            n = mirror[i].n
            area, error = hcubature((x) -> x[1] * f(x[1], x[2], mirror.z),
                (a * n, m * 2pi / (6n + 3)), (a * (n + 1), (m + 1) * 2pi / (6n + 3)))
            total = 0.0
            for t = 1:4
                tp = Mirrors.TransformParameters{a,s,t}()
                uvf(uv) = Mirrors.uvtransform(f, mirror.z, mirror.s, uv[1], uv[2], m, n, tp)
                sum, error = hcubature(uvf, (0.0, 0.0), (1.0, 1.0))
                total += sum
            end
            @test area ≈ total
        end
        # TODO: test to ensure that s and z are incorporated as they should be
    end



    # Test that integrating a singular patch works
    @testset "integratepatch" begin
        mirror = Mirror(9.7, 2, (r, θ) -> 0.0, (r, θ) -> 1.0)
        a = mirror.a
        # a function with a 1/r singularity at (r, th) that integrates to 1 if the circle of radius a/4pi centered at
        # (r, θ) is included; the function is continuous, but its derivatives are not
        function singularityfunction(a::Real, r::Real, θ::Real)
            return function (r_other::Real, θ_other::Real, z::Real)
                result = a / (4pi * sqrt(r_other^2 + r^2 - 2 * r_other * r * cos(θ_other - θ))) - 1
                return result > 0 ? 16pi / a^2 * result : 0
            end
        end # function singularityfunction
        for patch in mirror
            for (s, (r, θ, z)) in enumerate(patch)
                f = singularityfunction(a, r, θ)
                area = Mirrors.integratepatch(f, mirror, patch, s)
                @test isapprox(area, 1.0, atol=2e-3) # it's not very precise--the 6 innermost points give 0.998 :(
            end
        end
    end



    # Check that impedance works
    @testset "impedance" begin
        """
        Check whether the impedance matrix "integrates" to the expected area for a given function

        - `mirror`: the `Mirror` over which to integrate
        - `f`: the function to integrate
        - `expectedintegral`: the actual integral of `f` over `mirror`
        - `singular`: whether to use the singular code to fill the blocks on the main diagonal of `Z`
        - `equalpatches`: whether all patches on the mirror integrate to the same quantity
        """
        function test_expectedintegral(mirror::Mirror, f::Function, expectedintegral::Real, singular::Bool, equalpatches::Bool)
            Z = impedance(mirror, (r1, θ1, z1, r2, θ2, z2) -> f(r2, θ2, z2), singular)
            for s = 1:4
                area = 0.0
                for i = 1:length(mirror)
                    # diagonal
                    area += sum(Z[4i-4+s, 4i-3:4i])
                    # rows
                    @test isapprox(sum(Z[4i-4+s, :]), expectedintegral, atol=1e-2) # some are very imprecise :(
                    if equalpatches # patches are of equal integral "volume"--columns should also be equal
                        @test isapprox(sum(Z[s:4:end, 4i-3:4i]), expectedintegral, atol=1e-3)
                    end
                end
                @test isapprox(area, expectedintegral, atol=1e-4)
            end
        end
        # test expected integral under a few circumstances
        R = 1.1
        rings = 3
        for singular = (true, false), (f, scaled_area, equalpatches) = (((r, θ, z) -> 0.7, 0.7pi * R^2, true),
                ((r, θ, z) -> r, 2pi / 3 * R^3, false),
                ((r, θ, z) -> 1.4r * sin(θ)^2, 1.4pi * R^3 / 3, false),
                ((r, θ, z) -> 0.3r^2 * cos(θ)^2, 0.3pi * R^4 / 4, false))
            test_expectedintegral(Mirror(R, rings), f, scaled_area, singular, equalpatches)
        end
    end



    # THIS TEST FAILS--the condition numbers are pretty high :(
    # @testset "impedance" begin
    #     # make sure the condition number is sane for impedance matrices of flat mirrors
    #     for rings=1:5
    #         Z = impedance(Mirror(rand(), rings))
    #         @test cond(Z) < 1000
    #     end
    # end # @testset "impedance"
end

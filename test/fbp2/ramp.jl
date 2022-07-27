# test/fbp2/ramp.jl


using Sinograms: SinoPar, SinoMoj, SinoFanFlat, SinoFanArc, fbp_ramp
import Sinograms as SG
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred

include("../helper.jl")


@testset "ramp-unit" begin
    for ds_ in (2, 2f0, 2.0),
            dsd_ in (99, 99f0, 99.0),
            du in (1, 1mm)

        ds = ds_ * du
        dsd = dsd_ * du

        for n in 0:2
            f = @inferred SG._ramp_flat(n, ds)
            @test oneunit(f) == oneunit(1/ds^2)
            a = @inferred SG._ramp_arc(n, ds, dsd)
            @test oneunit(a) == oneunit(1/ds^2)
        end

        N = 6
        ha, na = @inferred SG.ramp_arc(N, ds, dsd)
        hf, nf = @inferred SG.ramp_flat(N, ds)

        n = 0:2
        fun = n -> SG._ramp_flat.(n, ds)
        fun(n)
        f = @inferred fun(n)
        fun = n -> SG._ramp_arc.(n, ds, dsd)
        a = @inferred fun(n)
    end
end


@testset "fbp-ramp" begin
    N = 8
    for shape in (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)
        isinteractive() && (@show shape)
        for d_ in (2, 2f0, 2.0),
                du in (1, 1mm),
                orbit_ in (90, 90f0, 90.0),
                ou in (1, 1°)

            sg = shape( ; d = d_*du, orbit=orbit_*ou)
            h, n = @inferred fbp_ramp(sg, N)
            @test h isa Vector
            @test length(h) == N
            @test length(n) == N
        end
    end
end

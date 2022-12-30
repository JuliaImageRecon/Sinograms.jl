#=
test/fbp2/parker.jl
=#

using Sinograms: SinoPar, SinoFanArc, SinoMoj
using Sinograms: CtFanArc
using Sinograms: angles, _ar, _gamma, _gamma_max, _orbit_short
using Sinograms: parker_weight, parker_weight_fan_short, parker_weight_par, dims
using Test: @test, @test_throws, @testset, @inferred


@testset "parker-par" begin
    rg = SinoPar()
    pw = @inferred parker_weight_par(rg.orbit, angles(rg))
    @test pw == ones(1,1)

    pw = @inferred parker_weight(rg)

    rg = SinoPar( ; orbit = 90)
    pw = @inferred parker_weight(rg) # warns

    rg = SinoPar( ; orbit = 361)
    @test_throws ErrorException parker_weight(rg)

    rg = SinoPar( ; orbit = 270)
    pw = @inferred parker_weight(rg)
    @test size(pw) == (1, rg.na)
    @test pw isa Matrix{Float32}
    @test pw[1] == 0 # todo: suboptimal?
    @test sum(pw)/length(pw) ≈ 1

    rg = SinoMoj()
    pw = @inferred parker_weight(rg)
    @test pw == ones(1,1)
end


@testset "parker-fan" begin
    rg = SinoFanArc(Val(:ge1), orbit=:short)
    pw = @inferred parker_weight_fan_short(
        rg.nb, rg.na, rg.orbit, _orbit_short(rg),
        _ar(rg), _gamma(rg), _gamma_max(rg),
    )
    @test pw isa Matrix{Float32}
    tmp = sum(pw; dims=2)
    tmp = 1 - minimum(tmp)/maximum(tmp)
    @test abs(tmp) < 0.002 # should be nearly uniform

    rg = SinoFanArc(; orbit = 360)
    pw = @inferred parker_weight(rg)
    @test pw isa Matrix{Float32}
    @test pw == ones(1,1)

    rg = SinoFanArc( ; orbit = 90)
    pw = @inferred parker_weight(rg) # warns
    @test pw isa Matrix{Float32}
end


@testset "parker-ct" begin
    cg = CtFanArc()
#   @inferred parker_weight_fan_short(cg)
    pw = @inferred parker_weight(cg)
    @test pw isa Array{Float32,3}
    @test size(pw) == (1,1,1)

    cg = CtFanArc(:short)
    pw = @inferred parker_weight(cg)
    @test pw isa Array{Float32,3}
    @test size(pw,2) == 1
end

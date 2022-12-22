#=
test/fbp2/parker.jl
=#

using Sinograms: SinoPar, SinoFanArc, SinoMoj
using Sinograms: parker_weight, parker_weight_fan, parker_weight_par, dims
using Test: @test, @test_throws, @testset, @inferred


@testset "parker-par" begin
    sg = SinoPar()
    pw = @inferred parker_weight_par(sg.orbit, sg.ad)
    @test pw == ones(sg.na)
    @test pw isa Vector{Float32}

    sg = SinoPar( ; orbit = 90)
    pw = @inferred parker_weight(sg) # warns

    sg = SinoPar( ; orbit = 361)
    @test_throws ErrorException parker_weight(sg)

    sg = SinoPar( ; orbit = 270)
    pw = @inferred parker_weight(sg)
    @test pw[1] == 0 # todo: suboptimal?
    @test sum(pw)/length(pw) ≈ 1

    sg = SinoMoj()
    pw = @inferred parker_weight(sg)
    @test pw isa Vector{Float32}
    @test pw[1] == 1
end


@testset "parker-fan" begin
    sg = SinoFanArc(Val(:ge1), orbit=:short)
    pw = @inferred parker_weight_fan(sg.nb, sg.na, sg.orbit, sg.orbit_short,
        sg.ar, sg.gamma, sg.gamma_max)
    @test pw isa Matrix{Float32}
    tmp = sum(pw; dims=2)
    tmp = 1 - minimum(tmp)/maximum(tmp)
    @test abs(tmp) < 0.002 # should be nearly uniform

    sg = SinoFanArc(; orbit=360)
    pw = @inferred parker_weight(sg)
    @test pw isa Matrix{Float32}
    @test pw == ones(dims(sg))

    sg = SinoFanArc( ; orbit = 90)
    pw = @inferred parker_weight(sg) # warns
    @test pw isa Matrix{Float32}
end

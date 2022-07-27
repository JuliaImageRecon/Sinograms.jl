#=
test/fbp2/parker.jl
=#

using Sinograms: SinoPar, SinoFanArc, SinoMoj
using Sinograms: parker_weight
using Test: @test, @test_throws, @testset, @inferred


@testset "parker" begin
    geoms = (SinoPar, SinoFanArc, SinoMoj)
    for geom in (geoms)
        sg = geom()
        pw = @inferred parker_weight(sg)
        @test pw isa Vector{Float32}
        @test pw[1] == 1
    end

    sg = SinoPar( ; orbit = 90)
    pw = @inferred parker_weight(sg)

    sg = SinoPar( ; orbit = 361)
    @test_throws String parker_weight(sg)

    sg = SinoPar( ; orbit = 270)
    pw = @inferred parker_weight(sg)
    @test pw[1] == 0 # todo: suboptimal?

    sg = SinoFanArc( ; orbit = 90)
    @test_throws String parker_weight(sg)
end

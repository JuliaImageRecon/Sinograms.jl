#=
test/fbp2/parker.jl
=#

using Sinograms: SinoPar, SinoFanArc, SinoMoj
using Sinograms: parker_weight, parker_weight_fan
using Test: @test, @test_throws, @testset, @inferred


@testset "parker" begin

    pw = @inferred parker_weight_fan(10, 360.)
    pw = @inferred parker_weight_fan(10, 99.)
    @test pw == zeros(10) # todo: until short-scan done

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
    @test sum(pw)/length(pw) ≈ 1

    sg = SinoFanArc( ; orbit = 90)
#   @test_throws String
    pw = @inferred parker_weight(sg)
    @test pw == zeros(sg.na) # todo: until short-scan done
end

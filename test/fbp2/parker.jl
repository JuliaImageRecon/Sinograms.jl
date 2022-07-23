#=
test/fbp2/parker.jl
=#

using Sinograms: SinoPar, SinoFanArc, SinoMoj
using Sinograms: parker_weight
using Test: @test, @testset, @inferred


@testset "parker" begin
    geoms = (SinoPar, SinoFanArc, SinoMoj)
    for geom in (geoms)
        sg = geom()
        pw = @inferred parker_weight(sg)
        @test pw isa Vector{Float32}
        @test pw[1] == 1
    end
end

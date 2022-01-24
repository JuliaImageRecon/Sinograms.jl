# runtests.jl

using Sinograms
using Test: @test, @testset, detect_ambiguities

include("helper.jl")

@testset "Sinograms" begin
    # todo: more

    include("fbp-par.jl")
    include("sino-geom.jl")

    @test length(detect_ambiguities(Sinograms)) == 0
end

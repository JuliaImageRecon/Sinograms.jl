# runtests.jl

using Sinograms
using Test: @test, @testset, detect_ambiguities

include("helper.jl")

@testset "Sinograms" begin
    # todo: more

    include("sino-geom.jl")
    include("fbp/window.jl")
    include("fbp2/ramp.jl")
    include("fbp-plan.jl")
    include("fbp-par.jl")
    include("zwart_powell.jl")

    @test length(detect_ambiguities(Sinograms)) == 0
end

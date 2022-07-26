# runtests.jl

using Sinograms
using Test: @test, @testset, detect_ambiguities

include("helper.jl")

include("fbp/window.jl")

include("fbp2/sino-geom.jl")
include("fbp2/parker.jl")
include("fbp2/ramp.jl")
include("fbp2/filter.jl")
include("fbp2/sino-plot.jl")
include("fbp2/back.jl")

include("fbp-plan.jl")
include("fbp-par.jl")
include("fbp-fan.jl")

include("sys2/zwart_powell.jl")

@testset "Sinograms" begin
    @test isempty(detect_ambiguities(Sinograms))
end

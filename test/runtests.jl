# runtests.jl

using Sinograms
using Test: @test, @testset, detect_ambiguities

include("helper.jl")

include("fbp/filter.jl")
include("fbp/ramp.jl")
include("fbp/window.jl")

include("fbp2/sino-geom.jl")
include("fbp2/parker.jl")
include("fbp2/sino-plot.jl")
include("fbp2/back2.jl")
include("fbp2/plan2.jl")

include("fbp-par.jl")
include("fbp-fan.jl")

include("fbp3/ct-geom.jl")
include("fbp3/plan3.jl")
include("fbp3/back3.jl")
include("fbp3/fdk.jl")
include("fbp3/ct-plot.jl")

include("sys2/zwart_powell.jl")

@testset "Sinograms" begin
    @test isempty(detect_ambiguities(Sinograms))
end

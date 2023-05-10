# runtests.jl

using Sinograms
using Test: @test, @testset, detect_ambiguities

include("units.jl")

include("geom/types.jl")
include("geom/sino-geom.jl")
include("geom/ct-geom.jl")

include("fbp/d-angle.jl")
include("fbp/filter.jl")
include("fbp/parker.jl")
include("fbp/ramp.jl")
include("fbp/window.jl")

include("fbp2/sino-plot.jl")
include("fbp2/back2.jl")
include("fbp2/plan2.jl")
include("fbp2/fbp-par.jl")
include("fbp2/fbp-fan.jl")

include("fbp3/plan3.jl")
include("fbp3/back3.jl")
include("fbp3/fdk.jl")
include("fbp3/ct-plot.jl")

include("sys2/zwart_powell.jl")
include("sys2/bdd_2d.jl")

@testset "Sinograms" begin
    @test isempty(detect_ambiguities(Sinograms))
end

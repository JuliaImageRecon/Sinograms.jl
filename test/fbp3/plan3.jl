#=
test/fbp3/plan3.jl
=#

using Sinograms: plan_fbp, Window, FDKplan
using Sinograms: CtFanArc, CtFanFlat, angles
using Sinograms: _fdk_weights, fdk_weight_cyl, _ar
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred

@testset "fbp3/plan3" begin
    rg = @inferred CtFanArc()
    ig = @inferred ImageGeom( ; dims=(32,30,10) )

    tmp = @inferred fdk_weight_cyl(rg)
    ar = @inferred _ar(rg)
    tmp = @inferred _fdk_weights(rg)

    plan = @inferred FDKplan(rg, ig, ones(rg.na), ones(2rg.ns))
    @test plan isa FDKplan

    plan = @inferred plan_fbp(rg, ig)
    @test plan isa FDKplan

    show(isinteractive() ? stdout : devnull, MIME("text/plain"), plan)

    rg = @inferred CtFanArc(:short)
    plan = @inferred plan_fbp(rg, ig)
    @test plan isa FDKplan
end

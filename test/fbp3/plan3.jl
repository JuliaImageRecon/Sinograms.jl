#=
test/fbp3/plan3.jl
=#

using Sinograms: plan_fbp, Window, FDKplan
using Sinograms: CtFanArc, CtFanFlat, angles
using Sinograms: _view_weights, _fdk_weights, fdk_weight_cyl
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred

#   @code_warntype angles(cg) # todo - needs work due to get property
#   @code_warntype fdk_weight_cyl(cg) # todo

@testset "fbp3/plan3" begin
    cg = CtFanArc()
    ig = ImageGeom( ; dims=(32,30,10) )

    tmp = @inferred fdk_weight_cyl(cg)
    ar = cg.ar
    da = @inferred _view_weights(ar)
    tmp = @inferred _fdk_weights(cg)

    plan = @inferred FDKplan(cg, ig, ones(cg.na), ones(2cg.ns))
    @test plan isa FDKplan

    plan = @inferred plan_fbp(cg, ig)
    @test plan isa FDKplan

    show(isinteractive() ? stdout : devnull, MIME("text/plain"), plan)

    cg = CtFanArc(:short)
    plan = @inferred plan_fbp(cg, ig)
    @test plan isa FDKplan
end

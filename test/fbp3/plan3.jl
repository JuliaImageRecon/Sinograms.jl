#=
test/fbp3/plan3.jl
=#

using Sinograms: plan_fbp, Window, FDKplan
using Sinograms: CtFanArc, CtFanFlat
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


@testset "fbp3/plan3" begin
    cg = CtFanArc()
    ig = ImageGeom( ; dims=(32,30,10) )

    plan = @inferred FDKplan(cg, ig, ones(cg.na), ones(2cg.ns))
    @test plan isa FDKplan

    plan = @inferred plan_fbp(cg, ig)
    @test plan isa FDKplan

    show(isinteractive() ? stdout : devnull, MIME("text/plain"), plan)
end

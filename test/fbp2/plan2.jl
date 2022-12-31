#=
test/fbp2/plan2.jl
=#

using Sinograms: plan_fbp, Window, FBPNormalPlan
using Sinograms: SinoPar, SinoFanArc
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


@testset "fbp2/plan2" begin
    rg = SinoPar()
    ig = ImageGeom()

    plan = @inferred FBPNormalPlan(rg, ig, ones(rg.na), ones(2rg.nb))
    @test plan isa FBPNormalPlan

    plan = @inferred plan_fbp(rg, ig)
    @test plan isa FBPNormalPlan

    show(isinteractive() ? stdout : devnull, MIME("text/plain"), plan)

    rg = @inferred SinoFanArc(:short)
    plan = @inferred plan_fbp(rg, ig)
    @test plan isa FBPNormalPlan
end

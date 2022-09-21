#=
test/fbp-plan.jl
=#

using Sinograms: plan_fbp, Window, FBPNormalPlan
using Sinograms: SinoPar
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


@testset "fbp-plan" begin
    sg = SinoPar()
    ig = ImageGeom()
    plan = @inferred FBPNormalPlan(sg, ig, ones(sg.na), ones(2sg.nb))
    @test plan isa FBPNormalPlan
    plan = @inferred plan_fbp(sg, ig)
    @test plan isa FBPNormalPlan
end

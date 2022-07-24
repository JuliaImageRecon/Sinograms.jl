#=
fbp-plan.jl
=#

using Sinograms: plan_fbp, Window, NormalPlan
using Sinograms: SinoPar
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


@testset "fbp-plan" begin
    sg = SinoPar()
    ig = ImageGeom()
    plan = @inferred NormalPlan(sg, ig, Window(), ones(sg.na))
    @test plan isa NormalPlan
    plan = @inferred plan_fbp(sg, ig)
    @test plan isa NormalPlan
end

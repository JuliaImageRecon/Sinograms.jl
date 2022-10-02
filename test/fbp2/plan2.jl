#=
test/fbp2/plan2.jl
=#

using Sinograms: plan_fbp, Window, FBPNormalPlan
using Sinograms: SinoPar
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


@testset "fbp2/plan2" begin
    sg = SinoPar()
    ig = ImageGeom()

    plan = @inferred FBPNormalPlan(sg, ig, ones(sg.na), ones(2sg.nb))
    @test plan isa FBPNormalPlan

    plan = @inferred plan_fbp(sg, ig)
    @test plan isa FBPNormalPlan

    show(isinteractive() ? stdout : devnull, MIME("text/plain"), plan)
end

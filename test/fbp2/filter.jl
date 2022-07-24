# test/fbp2/filter.jl


using Sinograms: SinoPar, fbp_filter
using Unitful: mm
using Test: @test, @testset, @test_throws, @inferred

include("../helper.jl")

@testset "fbp_filter" begin
    sg = SinoPar(; d = 2f0mm)
    Hk = @NOTinferred fbp_filter(sg) # todo: can't infer due to fbp_ramp?

    @test Hk isa Vector
    @test length(Hk) == 2 * sg.nb
end

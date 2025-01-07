# ext0.jl test extension stubs *before* loading the extension

using Sinograms: sino_plot_rays, sino_geom_plot!
using Sinograms: ct_geom_plot2!, ct_geom_plot3
using Test: @test_throws, @testset

@testset "ext0" begin
    @test_throws String sino_plot_rays()
    @test_throws String sino_geom_plot!()
    @test_throws String ct_geom_plot2!()
    @test_throws String ct_geom_plot3()
end


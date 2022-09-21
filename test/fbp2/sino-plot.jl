# test/fbp2/sino-plot.jl

# include("helper.jl")

using Plots: plot, gui, Plot
using Unitful: mm #, °
using UnitfulRecipes
#using MIRTjim: prompt
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
using Sinograms: sino_plot_rays, sino_geom_plot!
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


function _test_sino_plot_rays(geo ; nb::Int = 60, na::Int = 40)
    d, orbit = 2, 180
#   d, orbit = 2mm, 180.0°
    sg = geo(; d, orbit, nb, na)
    p1 = sino_plot_rays(sg)
    @test p1 isa Plot

    plot()
    p2 = sino_geom_plot!(sg)
    @test p2 isa Plot

    ig = ImageGeom()
    p3 = sino_geom_plot!(sg; ig)
    @test p2 isa Plot

    plot(p1, p2, p3)

    true
end


@testset "sino_plot" begin
#   sg = SinoFan(Val(:ge1))
    sg_list = (SinoPar, SinoFanArc, SinoFanFlat, SinoMoj)
    for geo in sg_list
        @test _test_sino_plot_rays(geo)
        gui()
#       prompt()
    end
end

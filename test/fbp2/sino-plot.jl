# test/fbp2/sino-plot.jl

# include("helper.jl")

using Plots: plot, gui, Plot
using Unitful: mm #, °
using UnitfulRecipes
#using MIRTjim: prompt
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
using Sinograms: sino_plot_rays, sino_geom_plot!
using Test: @test, @testset, @inferred


function _test_sino_plot_rays(geo ; nb::Int = 60, na::Int = 40)
    d, orbit = 2, 180
#   d, orbit = 2mm, 180.0°
    sg = geo(; d, orbit, nb, na)
    p = sino_plot_rays(sg)
    @test p isa Plot

    plot()
    p = sino_geom_plot!(sg)
    @test p isa Plot
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

# test/fbp2/sino-plot.jl

using Plots: plot, gui, Plot
using Unitful: mm #, °
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat, SinoFan
using Sinograms: sino_plot_rays, sino_geom_plot!
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


function _test_sino_plot_rays(geo ; nb::Int = 60, na::Int = 20, d=2mm)
    orbit = 180.0 #°
    arg = geo() isa SinoFan ? (; dod=120mm, dsd=300mm, d=2d) : (;)
    sg = geo(; d, orbit, nb, na, arg...)
    p1 = sino_plot_rays(sg)
    @test p1 isa Plot

    plot()
    ig = ImageGeom( ; dims = (64, 62), deltas = (1,1) .* d)
    p2 = sino_geom_plot!(sg, ig)
    @test p2 isa Plot

    return plot(p1, p2, layout = (2,1))
end


@testset "sino_plot" begin
#   sg = SinoFan(Val(:ge1))
    sg_list = (SinoPar, SinoFanArc, SinoFanFlat, SinoMoj)
    p = Array{Any}(undef, length(sg_list))
    for (i, geo) in enumerate(sg_list)
        p[i] = _test_sino_plot_rays(geo)
        @test p[i] isa Plot
    end
    plot(p..., layout = (1,4))
    gui()
end

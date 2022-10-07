# test/fbp3/ct-plot.jl

using Plots: plot, gui, Plot
using Unitful: mm #, °
using Sinograms: CtPar, CtFanArc, CtFanFlat, CtFan
using Sinograms: ct_geom_plot2!
using Sinograms: ct_geom_plot3
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


function _test_ct_plot(geo ;
    ns::Int = 60,
    nt::Int = 40,
    na::Int = 20,
    d = 2mm,
)

    orbit = 180.0 #°
    arg = geo() isa CtFan ? (; dod=120mm, dsd=300mm, ds=2d) : (; ds=d)
    sg = geo(; ns, nt, na, orbit, arg...)

    plot()
    ig = ImageGeom( ; dims = (64, 62, 30), deltas = (1,1,1) .* d)
    p = ct_geom_plot2!(sg, ig)

    p = ct_geom_plot3(sg, ig)

    return p
end


@testset "fbp3/ct-plot" begin
#   sg = SinoFan(Val(:ge1))
    list = (CtPar, CtFanArc, CtFanFlat)
   list = (CtFanArc, CtFanFlat)
    p = Array{Any}(undef, length(list))
    for (i, geo) in enumerate(list)
        p[i] = _test_ct_plot(geo)
        @test p[i] isa Plot
    end
    plot(p..., layout = (1,3))
    gui()
end

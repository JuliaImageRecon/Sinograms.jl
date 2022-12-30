# test/fbp3/ct-plot.jl

using Unitful: mm #, °
using Plots: plot, plot!, gui, Plot
using Sinograms: CtPar, CtFanArc, CtFanFlat, CtFan
using Sinograms: CtSourceCircle, CtSourceHelix
using Sinograms: ct_geom_plot2!
using Sinograms: ct_geom_plot3
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


function _test_ct_plot(geo ;
    ns::Int = 60,
    nt::Int = 40,
    na::Int = 20,
    d = 2mm,
    orbit = 180.0, #°
    src = CtSourceCircle(),
)

    arg = geo() isa CtFan ? (; dod=120mm, dsd=300mm, ds=2d) : (; ds=d)
    rg = geo(; ns, nt, na, orbit, src, arg...)

    plot()
    ig = ImageGeom( ; dims = (64, 62, 30), deltas = (1,1,1) .* d)
    p = ct_geom_plot2!(rg, ig)
    p = ct_geom_plot3(rg, ig)
    plot!(p; title = nameof(typeof(rg)))

    return p
end


@testset "fbp3/ct-plot" begin
#   rg = CtFanArc(Val(:ge1))
    list = (CtFanArc, CtFanFlat) # CtPar
    p = Array{Any}(undef, length(list))
    for (i, geo) in enumerate(list)
        p[i] = _test_ct_plot(geo)
        @test p[i] isa Plot
    end
    q = _test_ct_plot(CtFanArc; src=CtSourceHelix(;pitch=1.5,source_z0=0mm))
    plot!(q ; title="Helix")
    plot(p..., q)
    gui()
end

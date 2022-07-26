# test/fbp2/back.jl

using ImageGeoms: ImageGeom
using Sinograms: SinoPar, SinoFanArc, SinoFanFlat, fbp_back
import Sinograms as SG
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred
using LazyGrids: ndgrid

include("../helper.jl")

#=
    ig = @inferred ImageGeom()
    sg = SinoPar()
    sg = SinoFanArc()
    back = @inferred SG.fbp_back_par(sg.ones, sg.ar, sg.ds, sg.offset,
        sg.rfov, ndgrid(axes(ig)...)..., ig.mask, 1)
    back = @inferred fbp_back(sg, ig, sg.ones)
=#

@testset "back-par" begin
    ig = @inferred ImageGeom()
    for geom in (SinoPar, SinoFanArc, SinoFanFlat)
        sg = geom()
        back = @inferred fbp_back(sg, ig, sg.ones)
        @test back isa Matrix
    end
end

# test/fbp2/back.jl

using ImageGeoms: ImageGeom, MaskCircle
using LazyGrids: ndgrid
using Sinograms: SinoPar, SinoFanArc, SinoFanFlat, fbp_back
using Sinograms: fbp_back_par, fbp_back_par!, fbp_back_par_xy
using Unitful: mm
using Test: @test, @testset, @inferred

#=
# verify new pixel-driven vs old broadcasting/allocating way

using BenchmarkTools
using Random: seed!
using Sinograms: fbp_back_par_old

#   ig = @inferred ImageGeom(MaskCircle(), dims=(32,30), deltas=(1mm,1mm) )
    ig = @inferred ImageGeom(dims=(32,30), deltas=(1mm,1mm) )
#   sg = SinoPar( ; nb = 12, d=2mm) # intentionally small FOV
    sg = SinoPar( ; nb = 36, d=1mm)

    seed!(0)
    xc, yc = axes(ig)
    sino = randn(Float32, sg.dim)
    f1 = sino -> fbp_back_par_old(sino, sg.ar, sg.ds, sg.offset, # sg.rfov,
         ndgrid(xc, yc)..., ig.mask ; warned=true)
    f2 = sino -> fbp_back_par(sino, sg.ar, sg.ds, sg.offset,
         xc, yc, ig.mask)
    sang = sin.(sg.ar)
    cang = cos.(sg.ar)
    f3! = (image, sino) -> fbp_back_par!(image, sino, sang, cang,
         sg.ds, sg.offset, xc, yc, ig.mask)
    b1 = f1(sino)
    b2 = f2(sino)
    b3 = zeros(eltype(b1), size(b1))
    f3!(b3, sino)
    @test b1 ≈ b2 ≈ b3

#=
    x = xc[10] / sg.d
    y = yc[10] / sg.d
    f4 = sino -> fbp_back_par_xy(sino, sang, cang, sg.w,
          x, y ; T = Float32)
    f4(sino); # warm up

    @btime f4($sino);
=#

    @btime f1($sino);
    @btime f2($sino);
    @btime f3!($b3, $sino);
=#


@testset "fbp2/back" begin
    ig = @inferred ImageGeom( dims=(32,30), deltas=(1mm,1mm) )
    for geom in (SinoPar, SinoFanArc, SinoFanFlat)
        sg = geom( ; nb = 12, d=2mm) # intentionally small FOV
        back = @inferred fbp_back(sg, ig, sg.ones * 3mm; ia_skip = 2)
        @test back isa Matrix
    end
end

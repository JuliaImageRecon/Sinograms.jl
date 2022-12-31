# test/fbp2/back2.jl

using ImageGeoms: ImageGeom, MaskCircle
using LazyGrids: ndgrid
using Sinograms: SinoPar, SinoFanArc, SinoFanFlat, fbp_back, ones
using Unitful: mm
using Test: @test, @testset, @inferred

#=
# verify new pixel-driven vs old broadcasting/allocating way

using BenchmarkTools
using Random: seed!
using Sinograms: fbp_back_par, fbp_back_par!, fbp_back_par_xy
using Sinograms: fbp_back_fan, fbp_back_fan!, fbp_back_fan_xy
using Sinograms: fbp_back_par_old
using Sinograms: fbp_back_fan_old

#   ig = @inferred ImageGeom(MaskCircle(), dims=(32,30), deltas=(1mm,1mm) )
    ig = @inferred ImageGeom(dims=(32,30), deltas=(1mm,1mm) )
#   rg = SinoPar( ; nb = 12, d=2mm) # intentionally small FOV
    rg = SinoPar( ; nb = 36, d=1mm)
    rg = SinoFanArc( ; nb = 36, d=1mm)

    seed!(0)
    xc, yc = axes(ig)
    sino = randn(Float32, rg.dim)

  if rg isa SinoPar
    f1 = sino -> fbp_back_par_old(sino, rg.ar, rg.ds, rg.offset, # rg.rfov,
         ndgrid(xc, yc)..., ig.mask ; warned=true)
    f2 = sino -> fbp_back_par(sino, rg.ar, rg.ds, rg.offset,
         xc, yc, ig.mask)
    sang = sin.(rg.ar)
    cang = cos.(rg.ar)
    f3! = (image, sino) -> fbp_back_par!(image, sino, sang, cang,
         rg.ds, rg.offset, xc, yc, ig.mask)

  else # fan
    is_arc = rg isa SinoFanArc
    f1 = sino -> fbp_back_fan_old(sino, rg.ar,
         rg.dsd, rg.dso, rg.source_offset, is_arc,
         rg.ds, rg.offset, # rg.rfov,
         ndgrid(xc, yc)..., ig.mask ; warned=true)
    f2 = sino -> fbp_back_fan(sino, rg.ar,
         rg.dsd, rg.dso, rg.source_offset, is_arc,
         rg.ds, rg.offset,
         xc, yc, ig.mask)
    sβ = sin.(rg.ar)
    cβ = cos.(rg.ar)
    f3! = (image, sino) -> fbp_back_fan!(image, sino, sβ, cβ,
         rg.dsd, rg.dso, rg.source_offset, is_arc,
         rg.ds, rg.offset, xc, yc, ig.mask)
  end

    b1 = f1(sino)
    b2 = f2(sino)
    b3 = zeros(eltype(b1), size(b1))
    f3!(b3, sino)
    @test b1 ≈ b2 ≈ b3

#=
    x = xc[10] / rg.d
    y = yc[10] / rg.d
    f4 = sino -> fbp_back_par_xy(sino, sang, cang, rg.w,
          x, y ; T = Float32)
    f4(sino); # warm up

    @btime f4($sino);
=#

    @btime f1($sino);
    @btime f2($sino);
    @btime f3!($b3, $sino);
=#


@testset "fbp2/back2" begin
    ig = @inferred ImageGeom( dims=(32,30), deltas=(1mm,1mm) )
    for geom in (SinoPar, SinoFanArc, SinoFanFlat)
        rg = geom( ; nb = 12, d=2mm) # intentionally small FOV
        back = @inferred fbp_back(rg, ig, ones(rg) * 3mm; ia_skip = 2)
        @test back isa Matrix
    end
end

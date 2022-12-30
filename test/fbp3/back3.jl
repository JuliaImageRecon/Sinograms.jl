# test/fbp3/back3.jl

using ImageGeoms: ImageGeom, MaskCircle
using Sinograms: CtFanArc, CtFanFlat, cbct_back, dims
using Unitful: mm, kg
using Test: @test, @testset, @inferred

@testset "fbp3/back3" begin
    ig = @inferred ImageGeom( dims=(32,30,10), deltas=(1mm,1mm,2mm) )
    for geom in (CtFanArc, CtFanFlat)
        # intentionally small FOV:
        rg = @inferred geom( ; ns = 12, nt = 10, ds=2mm, na = 8)
        proj = fill(1kg, dims(rg))
        back = @inferred cbct_back(fill(1kg, dims(rg)), rg, ig ; ia_skip = 2)
        T = typeof(1f0 * oneunit(eltype(proj)))
        @test back isa Array{T, 3}
    end
end

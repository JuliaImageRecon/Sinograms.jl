# test/fbp3/back3.jl

using ImageGeoms: ImageGeom, MaskCircle
using LazyGrids: ndgrid
using Sinograms: CtFanArc, CtFanFlat, cbct_back, dims
using Unitful: mm
using Test: @test, @testset, @inferred

@testset "fbp3/back3" begin
    ig = @inferred ImageGeom( dims=(32,30,10), deltas=(1mm,1mm,2mm) )
    for geom in (CtFanArc, CtFanFlat)
        rg = @inferred geom( ; ns = 12, nt = 10, ds=2mm, na = 8) # intentionally small FOV
        back = @inferred cbct_back(ones(dims(rg)) * 3mm, rg, ig ; ia_skip = 2)
        @test back isa Array
    end
end

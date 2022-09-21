# test/fbp2/back.jl

using ImageGeoms: ImageGeom
using Sinograms: SinoPar, SinoFanArc, SinoFanFlat, fbp_back
using Unitful: mm
using Test: @test, @testset, @inferred

@testset "fbp2/back" begin
    ig = @inferred ImageGeom( dims=(32,30), deltas=(1mm,1mm) )
    for geom in (SinoPar, SinoFanArc, SinoFanFlat)
        sg = geom( ; nb = 12, d=2mm) # intentionally small FOV
        back = @inferred fbp_back(sg, ig, sg.ones * 3mm)
        @test back isa Matrix
    end
end

# test/geom/types.jl

using Sinograms # many
import Sinograms as SG
using Unitful: mm
using Test: @testset, @test, @inferred

@testset "geom/source" begin
    for traj in (CtSourceCircle, CtSourceHelix)
        src = @inferred traj()
        @test src isa CtSource
        show(isinteractive() ? stdout : devnull, MIME("text/plain"), src)
    end
    src = @inferred CtSourceCircle()
end

@testset "geom2" begin
    for d in (1f0, 1mm)
        for geo in (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)
            @test (@inferred geo( ; d)) isa SinoGeom
        end
        @test (@inferred SinoFanArc(Val(:ge1) ; unit = d)) isa SinoGeom
    end
end

@testset "geom3" begin
    for ds in (1f0, 1mm)
        for geo in (CtPar, CtFanArc, CtFanFlat) # CtMoj
            @test (@inferred geo( ; ds)) isa CtGeom
        end
        @test (@inferred CtFanArc(Val(:ge1) ; unit = ds)) isa CtGeom
        src = @inferred CtSourceHelix( ; pitch = 0.5)
        st = @inferred CtFanArc(Val(:ge1) ; unit = ds, src)
    end
end

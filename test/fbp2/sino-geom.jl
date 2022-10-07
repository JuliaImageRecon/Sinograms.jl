#=
test/fbp2/sino-geom.jl
=#

using Sinograms: RealU
using Sinograms: SinoGeom, SinoParallel, SinoFan
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
using Sinograms: ones, zeros, angles, rays, downsample, oversample # values
using Sinograms: sino_w, sino_s, footprint_size
import Sinograms as SG
using ImageGeoms: ImageGeom
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


function values(sg::S) where {S <: SinoGeom}
#   Iterators.map(p -> getproperty(sg, p), fieldnames(S))
    [getproperty(sg, p) for p in fieldnames(S)]
end


function _test(geo ; d=2mm, orbit=180.0°)
    sg = @inferred geo(; d, orbit)
#   args = 128, 100, d, orbit, 0*orbit, 0, 2d
    args = values(sg)
    @inferred geo{Int,Float64}(args...)
    sg = @inferred geo(args...)
    S = typeof(sg)
    sd = @inferred downsample(sg, 2)
    @test sd isa S
    so = @inferred oversample(sg, 2)
    @test so isa S
    true
end

@testset "construct" begin
    d, orbit = 2mm, 180.0°
#   d, orbit = 2, 180
    @inferred SinoPar(; d, orbit)
    @inferred SinoMoj(; d, orbit)
    @inferred SinoFanArc(; d, orbit)
    @inferred SinoFanFlat(; d, orbit)
    @inferred SinoFan(Val(:ge1))

    for geo in (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)
        @test _test(geo)
    end

    sg = @inferred SinoFan(Val(:ge1) ; orbit=:short)
    @test sg isa SinoFanArc
end


@testset "taufun" begin
    geoms = (
        (SinoPar(; d = 2, na = 5), (1:3)),
        (SinoPar(; d = 2, na = 5), (1:3) * 1.0),
        (SinoPar(; d = 2mm, na = 5), (1:3)mm),
        (SinoPar(; d = 2.0mm, na = 5), (1:3)mm),
    )
    for tmp in geoms
        sg = tmp[1]
        x = tmp[2]
        @inferred SG.sino_geom_tau(sg, x, 2x)
        @inferred sg.taufun(x, 2*x)
    end
end


function _test_prop(sg; d = 2mm, orbit = 180.0°)
    sg.ad[2]
    @test sg.rfov isa RealU
    @test (@inferred sg.down(2)) isa SinoGeom
    @test (@inferred sg.over(2)) isa SinoGeom
    @test sg.dim isa Dims
    @test sg.w isa RealU
    @test sg.ones isa Array
    @test sg.zeros isa Array
    @test sg.dr isa RealU
    @test sg.ds isa RealU
    @test sg.r isa AbstractVector
    @test sg.s isa AbstractVector
    @test sg.ad isa AbstractVector
    @test sg.ar isa AbstractVector
    @test sg.xds isa AbstractVector
    @test sg.yds isa AbstractVector

    show(isinteractive() ? stdout : devnull, sg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), sg)
    @test (@inferred ones(sg)) isa AbstractMatrix{Float32}
    @test (@inferred zeros(sg)) isa AbstractMatrix{Float32}
    @test (@inferred sino_w(sg)) isa Number
    @test (@inferred sino_s(sg)) isa AbstractVector
    @test (@inferred angles(sg)) isa AbstractVector
    rs = @inferred rays(sg)
    @test rs isa Tuple
    @test (@inferred propertynames(sg)) isa NTuple

    γ = @inferred SG.sino_geom_gamma(sg)
    @test γ isa Union{Nothing, AbstractVector}

    if sg isa SinoFan
        @test sg.gamma isa AbstractVector
        @test sg.gamma_max isa RealU
        @test sg.orbit_short isa RealU
        @test sg.dsd isa RealU
        @test sg.dfs isa RealU
        @test sg.dso isa RealU
        @test sg.dod isa RealU
    end

    if sg isa SinoMoj
        @test sg.d_moj(0) isa RealU
        @test sg.d_ang isa AbstractVector # angular dependent d for :moj
    end

    @test sg.shape(vec(sg.ones)) == sg.ones
    x = (1:4) * oneunit(d)
    @inferred sg.taufun(x, 2*x)
    @test sg.unitv() isa Array
    @test sg.unitv(ib=1, ia=2) isa Array

    ig = ImageGeom( (5,7), (1,1) .* d)
    @test (@inferred footprint_size(sg, ig)) isa Float32

    true
end


sg_list = (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)

@testset "orbit-start" begin
    down = 8
    d, orbit = 2mm, 180.0°
    orbit_start = 20°
    for geo in sg_list
        sg = @inferred geo( ; d, orbit, orbit_start, down)
        @test _test_prop(sg)
    end
end

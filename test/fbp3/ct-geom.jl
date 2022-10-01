#=
test/fbp3/ct-geom.jl
=#

using Sinograms: RealU
using Sinograms: CtGeom, CtParallel, CtFan
using Sinograms: CtPar, CtFanArc, CtFanFlat
using Sinograms: dims, ones, zeros, angles, rays, downsample, oversample # values
import Sinograms as SG
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


function values(cg::S) where {S <: CtGeom}
#   Iterators.map(p -> getproperty(cg, p), fieldnames(S))
    [getproperty(cg, p) for p in fieldnames(S)]
end


function _test(geo ; ds=2mm, orbit=180.0°)
    cg = @inferred geo(; ds, orbit)
#   args = 128, 100, ds, orbit, 0*orbit, 0, 2ds
    args = values(cg)
#   @inferred geo{Int,Float64}(args...)
    cg = @inferred geo(args...)
    S = typeof(cg)
    cd = @inferred downsample(cg, 2)
    @test cd isa S
#   co = @inferred oversample(cg, 2)
#   @test co isa S
    true
end

@testset "construct" begin
    ds, orbit = 2mm, 180.0°
#   ds, orbit = 2, 180
    @inferred CtPar(; ds, orbit)
    @inferred CtFanArc(; ds, orbit)
    @inferred CtFanFlat(; ds, orbit)
    @test (@inferred SG.ct_geom_ge1()) isa CtFanArc

    for geo in (CtPar, CtFanArc, CtFanFlat)
        @test _test(geo)
    end

#   cg = @inferred CtFan(Val(:ge1) ; orbit=:short)
#   @test cg isa CtFanArc
end


#=
@testset "taufun" begin
    geoms = (
        (CtPar(; d = 2, na = 5), (1:3)),
        (CtPar(; d = 2, na = 5), (1:3) * 1.0),
        (CtPar(; d = 2mm, na = 5), (1:3)mm),
        (CtPar(; d = 2.0mm, na = 5), (1:3)mm),
    )
    for tmp in geoms
        cg = tmp[1]
        x = tmp[2]
        @inferred SG.ct_geom_taufun(cg, x, 2x)
        @inferred cg.taufun(x, 2*x)
    end
end
=#


function _test_prop(cg; d = 2mm, orbit = 180.0°)
    @test cg.ad isa AbstractVector
#   cg.rfov
#   @inferred cg.down(2)
    @test (@inferred downsample(cg,2)) isa CtGeom
#   @inferred cg.over(2)
#   cg.dims
    @test cg.ws isa Real
    @test cg.wt isa Real
#   cg.ones
#   cg.zeros
    @test cg.ds isa RealU
    @test cg.dt isa RealU
    @test cg.s isa AbstractVector
    @test cg.t isa AbstractVector
    @test cg.ad isa AbstractVector
    @test cg.ar isa AbstractVector
    @test cg.xds isa AbstractVector
    @test cg.yds isa AbstractVector
    @test cg.cone_angle isa Real
    @test cg.rfov isa RealU
    @test cg.zfov isa RealU
    @test cg.unitv() isa Array
    @test cg.unitv( ; is=1) isa Array

    show(isinteractive() ? stdout : devnull, cg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), cg)

    # common methods
    @test dims(cg) isa Dims{3}
    @test (@inferred angles(cg)) isa AbstractVector

    @test (@inferred SG.ct_geom_ws(cg)) isa SG.Toffset
    @test (@inferred SG.ct_geom_wt(cg)) isa SG.Toffset
    @test (@inferred SG.ct_geom_s(cg)) isa LinRange
    @test (@inferred SG.ct_geom_t(cg)) isa LinRange

#   @test (@inferred propertynames(cg)) isa NTuple # todo
    @test (propertynames(cg)) isa NTuple

    @test (@inferred ones(cg)) isa Array{Float32,3}
    @test (@inferred zeros(cg)) isa Array{Float32,3}

    if cg isa CtPar
        rs = @inferred rays(cg)
        @test rs isa Base.Iterators.ProductIterator
    else
        rs = rays(cg) # @NOTinferred todo
        @test rs isa Base.Generator
    end
	@test collect(rs) isa Array{<:Tuple}

#   γ = @inferred SG.ct_geom_gamma(cg)
#   @test γ isa Union{Nothing, AbstractVector}

#   if cg isa CtParallel
#       @test isnothing(cg.gamma)
#   end

    if cg isa CtFan
        @test cg.gamma isa AbstractVector
        @test cg.gamma_max isa Real
        @test cg.orbit_short isa Real
        @test cg.dsd isa RealU
        @test cg.dfs isa RealU
        @test cg.dso isa RealU
        @test cg.dod isa RealU
        @test cg.source_dz_per_view isa RealU
        @test cg.source_zs isa AbstractVector
        @test cg.gamma_max_abs isa Real
    end

    @test cg.shape(vec(ones(cg))) == ones(cg)

#=
    x = (1:4) * oneunit(d)
    @inferred cg.taufun(x, 2*x)
    cg.unitv()
    cg.unitv(ib=1, ia=2)
=#
    true
end


cg_list = (CtPar, CtFanArc, CtFanFlat)

@testset "orbit-start" begin
    down = 8
#   ds, orbit = 2mm, 180.0° # todo: degrees later
    ds, orbit = 2mm, 180f0
    orbit_start = 20 * oneunit(orbit)
    for geo in cg_list
        cg = @inferred geo( ; ds, orbit, orbit_start, down)
        @test _test_prop(cg)
    end
end

@testset "pitch" begin
    cg = CtFanArc( ; pitch = 2)
    @test cg.source_dz_per_view isa RealU
end

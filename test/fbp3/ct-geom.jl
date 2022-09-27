#=
test/fbp3/ct-geom.jl
=#

#using Sinograms: RealU
using Sinograms: CtGeom, CtParallel, CtFan
using Sinograms: CtPar, CtFanArc, CtFanFlat
using Sinograms: ones, zeros, angles, rays, downsample, oversample # values
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
#   @inferred CtFan(Val(:ge1))

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
    cg.ad[2]
#   cg.rfov
#   @inferred cg.down(2)
    @inferred downsample(cg,2)
#   @inferred cg.over(2)
#   cg.dims
    cg.ws
    cg.wt
#   cg.ones
#   cg.zeros
    cg.ds
    cg.dt
    cg.s
    cg.t
    cg.ad
    cg.ar
    cg.xds
    cg.yds

#=
#todo after debugging
    show(isinteractive() ? stdout : devnull, cg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), cg)
=#

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

#   γ = @inferred SG.ct_geom_gamma(cg)
#   @test γ isa Union{Nothing, AbstractVector}

#   if cg isa CtParallel
#       @test isnothing(cg.gamma)
#   end

    if cg isa CtFan
        @test cg.gamma isa AbstractVector
        cg.gamma_max
        @test cg.orbit_short isa Real
        cg.dsd
        cg.dfs
        cg.dso
#       cg.dod
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

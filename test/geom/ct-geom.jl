#=
test/geom/ct-geom.jl
=#

using Sinograms: RealU
using Sinograms: CtGeom, CtParallel, CtFan
using Sinograms: CtPar, CtFanArc, CtFanFlat
using Sinograms: dims, ones, zeros, angles, rays, axes, downsample, oversample
using Sinograms: _ds, _ar, _xds, _yds, _rfov, footprint_size, _orbit_short
using Sinograms: _s, _t, _ws, _wt
using Sinograms: _tau, _shape, _unitv
using Sinograms: _gamma, _gamma_max, _dfs, _dso
using Sinograms: _source_dz_per_view, _source_zs, _zfov, _cone_angle, CtSourceHelix
using ImageGeoms: ImageGeom
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


list = (CtPar, CtFanArc, CtFanFlat)
ge1 = (; kwargs...) -> CtFanArc(Val(:ge1) ; kwargs...)

function values(rg::R) where {R <: CtGeom}
#   Iterators.map(p -> getfield(rg, p), fieldnames(R))
    [getfield(rg, p) for p in fieldnames(R)]
end


function _test_construct3(geo ; ds=2mm, orbit=180.0°)
    rg = @inferred geo(; ds, orbit)
    args = values(rg)
    Td = typeof(ds)
    To = typeof(orbit)
    Ts = eltype(rg.src)

    @inferred geo{Td,To,Ts}(args...)

    rg = @inferred geo(args...)
    R = typeof(rg)
    sd = @inferred downsample(rg, 2)
    @test sd isa R
    so = @inferred oversample(rg, 2)
    @test so isa R
    true
end


@testset "construct" begin
    ds, orbit = 2mm, 180.0° # stress
    for geo in list
        @test (@inferred geo(; ds, orbit)) isa CtGeom
        @test _test_construct3(geo)
    end

    rg = @inferred ge1( ; unit = oneunit(ds), orbit=:short)
    @test rg isa CtFanArc

    rg = @inferred CtFanArc(:short)
    @test rg isa CtFanArc
    rg = @inferred CtFanFlat(:short)
    @test rg isa CtFanFlat
end


@testset "tau" begin
    rg = CtPar()
    x = 1:3
    @inferred _tau(rg, x, 2x)
    @inferred _tau(rg)(x, 2x)
    x = ones(1,1) # Array
    @inferred _tau(rg)(x, 2x)
end


function _test_prop(rg::CtGeom{Td,To}) where {Td,To}

    @test (@inferred _ws(rg)) isa Real
    @test (@inferred _wt(rg)) isa Real

    @test (@inferred _s(rg)) isa AbstractVector
    @test (@inferred _t(rg)) isa AbstractVector

    @test (@inferred _ar(rg)) isa AbstractVector
    @test (@inferred _xds(rg)) isa AbstractVector
    @test (@inferred _yds(rg)) isa AbstractVector

    @test (@inferred _rfov(rg)) isa RealU
    @test (@inferred _zfov(rg)) isa RealU

    @test (@inferred _unitv(rg)) isa Array
    @test (@inferred _unitv(rg, (1,2,3))) isa Array

    @test (@inferred _source_dz_per_view(rg)) isa RealU
    @test (@inferred _source_zs(rg)) isa AbstractVector

    show(isinteractive() ? stdout : devnull, rg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), rg)

    # common methods
    D = length(dims(rg))
    @test (@inferred dims(rg)) isa Dims{D}
    @test (@inferred axes(rg)) isa Tuple
    @test (@inferred angles(rg)) isa AbstractVector
    @test (@inferred ones(rg)) isa Array{Float32,D}
    @test (@inferred zeros(rg)) isa Array{Float32,D}
    @test (@inferred _shape(rg, vec(ones(rg)))) == ones(rg)
    @test (@inferred _shape(rg, vec(ones(dims(rg)...,2)), :)) == ones(dims(rg)...,2)

    if rg isa CtFan
        @test (@inferred _ds(rg)) isa RealU

        @test rg.dsd isa RealU
        @test rg.dod isa RealU
        @test (@inferred _dfs(rg)) isa RealU
        @test (@inferred _dso(rg)) isa RealU

        @test (@inferred _orbit_short(rg)) isa RealU
        @test (@inferred _cone_angle(rg)) isa RealU

        @test (@inferred _gamma(rg)) isa AbstractVector
        @test (@inferred _gamma(rg, _s(rg)[1:2:end])) isa AbstractVector
        @test (@inferred _gamma_max(rg)) isa Real
    end


    if rg isa CtPar
        rs = @inferred rays(rg)
        @test rs isa Base.Iterators.ProductIterator
    else
        rs = @inferred rays(rg)
        @test rs isa Base.Generator{<:Base.Iterators.ProductIterator}
    end

    @test collect(rs) isa Array{<:Tuple}

    x = (1:4) * oneunit(Td)
    @inferred _tau(rg)(x, 2x)

    ig = ImageGeom((5,5,5), (1,1,1) .* oneunit(Td))
    @test (@inferred footprint_size(rg, ig)) isa Float32

    true
end


@testset "orbit-start" begin
    ds, orbit = 2mm, 180.0° # stress
    orbit_start = 20 * oneunit(orbit)
    for geo in list
        rg = @inferred geo( ; ns = 6, nt = 8, na = 5, ds, orbit, orbit_start)
        @test _test_prop(rg)
    end
end

@testset "pitch" begin
    src = @inferred CtSourceHelix(; pitch = 2)
    rg = @inferred CtFanArc( ; src)
    @test _source_dz_per_view(rg) isa RealU
end

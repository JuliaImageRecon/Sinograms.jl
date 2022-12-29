#=
test/geom/sino-geom.jl
=#

using Sinograms: RealU
using Sinograms: SinoGeom, SinoParallel, SinoFan
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
using Sinograms: dims, ones, zeros, angles, rays, axes, downsample, oversample
using Sinograms: _ar, _xds, _yds, _rfov, footprint_size, _orbit_short
using Sinograms: _r, _s, _w, _tau, _taufun, _d_moj, _d_ang
# using Sinograms: _ds, _dr
using Sinograms: _shape, _unitv
using Sinograms: _gamma, _gamma_max, _gamma_max_abs, _dfs, _dso
using ImageGeoms: ImageGeom
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


list = (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)
ge1 = (; kwargs...) -> SinoFanArc(Val(:ge1) ; kwargs...)

function values(rg::R) where {R <: SinoGeom}
#   Iterators.map(p -> getfield(rg, p), fieldnames(R))
    [getfield(rg, p) for p in fieldnames(R)]
end


function _test_construct(geo ; d=2mm, orbit=180.0°)
    rg = @inferred geo(; d, orbit)
    args = values(rg)
    Td = eltype(d)
    To = eltype(orbit)

    @inferred geo{Td,To}(args...)

    rg = @inferred geo(args...)
    R = typeof(rg)
    sd = @inferred downsample(rg, 2)
    @test sd isa R
    so = @inferred oversample(rg, 2)
    @test so isa R
    true
end


@testset "construct" begin
    d, orbit = 2mm, 180.0° # stress
    for geo in list
        @test (@inferred geo(; d, orbit)) isa SinoGeom
        @test _test_construct(geo)
    end

    rg = @inferred ge1( ; unit = oneunit(d), orbit=:short)
    @test rg isa SinoFanArc

    rg = @inferred SinoFanArc(:short)
    @test rg isa SinoFanArc
    rg = @inferred SinoFanFlat(:short)
    @test rg isa SinoFanFlat
end


@testset "taufun" begin
    geoms = (
        (SinoPar(; d = 2, na = 5), (1:3)),
        (SinoPar(; d = 2, na = 5), (1:3) * 1.0),
        (SinoPar(; d = 2mm, na = 5), (1:3)mm),
        (SinoPar(; d = 2.0mm, na = 5), (1:3)mm),
    )
    for (rg, x) in geoms
        @inferred _tau(rg, x, 2x)
        @inferred _taufun(rg)(x, 2x)
    end

    rg = SinoPar()
    x = ones(1,1) # Array
    @inferred _taufun(rg)(x, 2x)
end


function _test_prop(
    rg::SinoGeom{Td,To} ;
    d = 2*oneunit(Td),
    orbit = 180.0*oneunit(Td),
) where {Td,To}

    @test (@inferred _w(rg)) isa RealU
#   @test (@inferred _dr(rg)) isa RealU
#   @test (@inferred _ds(rg)) isa RealU
    @test (@inferred _r(rg)) isa AbstractVector
    @test (@inferred _s(rg)) isa AbstractVector

    @test (@inferred _ar(rg)) isa AbstractVector
    @test (@inferred _xds(rg)) isa AbstractVector
    @test (@inferred _yds(rg)) isa AbstractVector

    @test (@inferred _rfov(rg)) isa RealU

    @test (@inferred _unitv(rg)) isa Array
    @test (@inferred _unitv(rg, (1,2))) isa Array

    show(isinteractive() ? stdout : devnull, rg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), rg)

    # common methods
    D = length(dims(rg))
    @test (@inferred dims(rg)) isa Dims{D}
    @test (@inferred axes(rg)) isa Tuple
    @test (@inferred angles(rg)) isa AbstractVector
    @test (@inferred ones(rg)) isa Array{Float32,D}
    @test (@inferred zeros(rg)) isa Array{Float32,D}
    @test _shape(rg, vec(ones(rg))) == ones(rg) # NOTinferred

    if rg isa SinoFan
        @test rg.dsd isa RealU
        @test rg.dod isa RealU
        @test (@inferred _dfs(rg)) isa RealU
        @test (@inferred _dso(rg)) isa RealU

        @test (@inferred _orbit_short(rg)) isa RealU

        @test (@inferred _gamma(rg)) isa AbstractVector
        @test (@inferred _gamma(rg, _s(rg)[1:2:end])) isa AbstractVector
        @test (@inferred _gamma_max(rg)) isa Real
        @test (@inferred _gamma_max_abs(rg)) isa Real
    end

    if rg isa SinoMoj
        @test _d_moj(rg)(0) isa RealU
        @test _d_ang(rg) isa AbstractVector # angular dependent d for :moj
    end

    if rg isa SinoPar
        rs = @inferred rays(rg)
        @test rs isa Iterators.ProductIterator
    elseif rg isa SinoMoj
        rs = @inferred rays(rg)
        @test rs isa Base.Generator{<:Base.Iterators.ProductIterator}
    else
        rs = rays(rg) # @NOTinferred because of "fun" closure
        @test rs isa Base.Generator{<:Base.Iterators.ProductIterator}
    end

    @test collect(rs) isa Array{<:Tuple} # @NOTinferred

    ig = ImageGeom( (5,7), (1,1) .* d)
    @test (@inferred footprint_size(rg, ig)) isa Float32

    true
end


@testset "orbit-start" begin
    d, orbit = 2mm, 180.0° # stress
    orbit_start = 20 * oneunit(orbit)
    for geo in list
        rg = @inferred geo( ; nb = 6, na = 5, d, orbit, orbit_start)
        @test _test_prop(rg)
    end
end

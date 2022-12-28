#=
test/geom/ct-geom.jl
=#

using Sinograms: RealU, CtSourceHelix
using Sinograms: CtGeom, CtParallel, CtFan
using Sinograms: CtPar, CtFanArc, CtFanFlat
using Sinograms: dims, ones, zeros, angles, rays, axes, downsample, oversample
using Sinograms: _s, _t, _ws, _wt, _ar, _xds, _yds, _rfov, _zfov
using Sinograms: _shape, _unitv
using Sinograms: _gamma, _gamma_max, _gamma_max_abs, _dfs, _dso
using Sinograms: footprint_size, _source_dz_per_view, _source_zs
using Sinograms: _orbit_short, _cone_angle
import Sinograms as SG
using ImageGeoms: ImageGeom
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


list = (CtPar, CtFanArc, CtFanFlat)
ge1 = (; kwargs...) -> CtFanArc(Val(:ge1) ; kwargs...)

function values(st::S) where {S <: CtGeom}
#   Iterators.map(p -> getproperty(st, p), fieldnames(S))
    [getproperty(st, p) for p in fieldnames(S)]
end


function _test_construct(geo ; ds=2mm, orbit=180.0°)
    st = @inferred geo(; ds, orbit)
    args = values(st)
    Td = eltype(ds)
    To = eltype(orbit)
    Ts = eltype(st.src)
    @inferred geo{Td,To,Ts}(args...)
    st = @inferred geo(args...)
    S = typeof(st)
    sd = @inferred downsample(st, 2)
    @test sd isa S
    so = @inferred oversample(st, 2)
    @test so isa S
    true
end

@testset "construct" begin
    ds, orbit = 2mm, 180.0° # stress
    for geo in list
        @test (@inferred geo(; ds, orbit)) isa CtGeom
        @test _test_construct(geo)
    end

    st = @inferred ge1( ; unit = oneunit(ds), orbit=:short)
    @test st isa CtFanArc

    st = CtFanArc(:short) # @NOTinferred
    @test st isa CtFanArc
    st = CtFanFlat(:short) # @NOTinferred
    @test st isa CtFanFlat
end


function _test_prop(
    st::CtGeom{Td,To} ;
    d = 2*oneunit(Td),
    orbit = 180.0*oneunit(Td),
) where {Td,To}

    @test st.ds isa RealU
    @test st.dt isa RealU
#=
    @test st.ws isa Real
    @test st.wt isa Real
    @test st.s isa AbstractVector
    @test st.t isa AbstractVector
=#

    @test (@inferred _ws(st)) isa Real
    @test (@inferred _wt(st)) isa Real
    @test (@inferred _s(st)) isa AbstractVector
    @test (@inferred _t(st)) isa AbstractVector

    @test (@inferred _ar(st)) isa AbstractVector
    @test (@inferred _xds(st)) isa AbstractVector
    @test (@inferred _yds(st)) isa AbstractVector

    @test (@inferred _rfov(st)) isa RealU
    @test (@inferred _zfov(st)) isa RealU

    @test (@inferred _unitv(st)) isa Array
    @test (@inferred _unitv(st, (1,2,3))) isa Array

    @test (@inferred _source_dz_per_view(st)) isa RealU
    @test (@inferred _source_zs(st)) isa AbstractVector

    show(isinteractive() ? stdout : devnull, st)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), st)

#= todo xx
    @test st.ad isa AbstractVector
    @test st.ar isa AbstractVector
    @test st.xds isa AbstractVector
    @test st.yds isa AbstractVector
    @test st.rfov isa RealU
    @test st.zfov isa RealU
    @test st.unitv() isa Array
    @test st.unitv((1,2,3)) isa Array
=#

    # common methods
    D = length(dims(st))
    @test (@inferred dims(st)) isa Dims{D}
    @test (@inferred axes(st)) isa Tuple
    @test (@inferred angles(st)) isa AbstractVector
    @test (@inferred ones(st)) isa Array{Float32,D}
    @test (@inferred zeros(st)) isa Array{Float32,D}
    @test _shape(st, vec(ones(st))) == ones(st) # NOTinferred

#= todo xx
    @test (@inferred SG.ct_geom_ws(st)) isa SG.Toffset
    @test (@inferred SG.ct_geom_wt(st)) isa SG.Toffset
    @test (@inferred SG.ct_geom_s(st)) isa LinRange
    @test (@inferred SG.ct_geom_t(st)) isa LinRange

    @test st.source_dz_per_view isa RealU
    @test st.source_zs isa AbstractVector

    @test (@inferred propertynames(st)) isa NTuple
=#

    if st isa CtFan
        @test st.dsd isa RealU
        @test st.dod isa RealU
        @test (@inferred _dfs(st)) isa RealU
        @test (@inferred _dso(st)) isa RealU

        @test (@inferred _orbit_short(st)) isa RealU
        @test (@inferred _cone_angle(st)) isa RealU

        @test (@inferred _gamma(st)) isa AbstractVector
        @test (@inferred _gamma(st, _s(st)[1:2:end])) isa AbstractVector
        @test (@inferred _gamma_max(st)) isa Real
        @test (@inferred _gamma_max_abs(st)) isa Real
#= todo xx
        @test st.gamma isa AbstractVector
        @test st.gamma_s(st.s) isa AbstractVector
        @test st.gamma_max isa RealU
        @test st.gamma_max_abs isa Real
        @test st.orbit_short isa RealU
        @test st.dso isa RealU
        @test st.dfs isa RealU

        @test st.cone_angle isa Real
=#
    end


    if st isa CtPar
        rs = @inferred rays(st)
        @test rs isa Base.Iterators.ProductIterator
    else
        rs = @inferred rays(st) # @NOTinferred because of "fun" closure
        @test rs isa Base.Generator{<:Base.Iterators.ProductIterator}
    end

    @test collect(rs) isa Array{<:Tuple} # @NOTinferred


#=
    x = (1:4) * oneunit(d)
    @inferred st.taufun(x, 2*x)
=#

    ig = ImageGeom((9,9,9), (1,1,1) .* st.ds)
    @test (@inferred footprint_size(st, ig)) isa Float32

    true
end


@testset "orbit-start" begin
    ds, orbit = 2mm, 180.0° # stress
#   ds, orbit = 2mm, 180f0
    orbit_start = 20 * oneunit(orbit)
    for geo in list
        st = @inferred geo( ; ns = 6, nt = 8, na = 5, ds, orbit, orbit_start)
        @test _test_prop(st)
    end
end

@testset "pitch" begin
    src = @inferred CtSourceHelix(; pitch = 2)
    st = @inferred CtFanArc( ; src)
    @test _source_dz_per_view(st) isa RealU
end

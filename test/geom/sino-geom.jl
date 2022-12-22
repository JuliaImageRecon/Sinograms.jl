#=
test/geom/sino-geom.jl
=#

using Sinograms: RealU
using Sinograms: SinoGeom, SinoParallel, SinoFan
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
using Sinograms: dims, ones, zeros, angles, rays, downsample, oversample
using Sinograms: sino_w, sino_s, footprint_size
import Sinograms as SG
using ImageGeoms: ImageGeom
using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


list = (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)
ge1 = (; kwargs...) -> SinoFanArc(Val(:ge1) ; kwargs...)

function values(st::S) where {S <: SinoGeom}
#   Iterators.map(p -> getproperty(st, p), fieldnames(S))
    [getproperty(st, p) for p in fieldnames(S)]
end


function _test_construct(geo ; d=2mm, orbit=180.0°)
    st = @inferred geo(; d, orbit)
    args = values(st)
    Td = eltype(d)
    To = eltype(orbit)
    @inferred geo{Td,To}(args...)
    st = @inferred geo(args...)
    S = typeof(st)
    sd = @inferred downsample(st, 2)
    @test sd isa S
    so = @inferred oversample(st, 2)
    @test so isa S
    true
end

@testset "construct" begin
    d, orbit = 2mm, 180.0° # stress
    for geo in list
        @test (@inferred geo(; d, orbit)) isa SinoGeom
        @test _test_construct(geo)
    end

    st = @inferred ge1( ; unit = oneunit(d), orbit=:short)
    @test st isa SinoFanArc

    st = SinoFanArc(:short) # @NOTinferred
    @test st isa SinoFanArc
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


function _test_prop(
    st::SinoGeom{Td,To} ;
    d = 2*oneunit(Td),
    orbit = 180.0*oneunit(Td),
) where {Td,To}

    @test st.w isa RealU
    @test st.dr isa RealU
    @test st.ds isa RealU
    @test st.r isa AbstractVector
    @test st.s isa AbstractVector

    @test st.ad isa AbstractVector
    @test st.ar isa AbstractVector
    @test st.xds isa AbstractVector
    @test st.yds isa AbstractVector
    @test st.rfov isa RealU
    @test st.unitv() isa Array
    @test st.unitv((1,2)) isa Array

    show(isinteractive() ? stdout : devnull, st)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), st)

    # common methods
    D = length(dims(st))
    @test (@inferred dims(st)) isa Dims{D}
    @test (@inferred angles(st)) isa AbstractVector
    @test (@inferred ones(st)) isa Array{Float32,D}
    @test (@inferred zeros(st)) isa Array{Float32,D}
    @test st.shape(vec(ones(st))) == ones(st)

    @test (@inferred sino_w(st)) isa SG.Toffset
    @test (@inferred sino_s(st)) isa AbstractVector

    @test (@inferred propertynames(st)) isa NTuple

    if st isa SinoFan
        γ = @inferred SG._geom_gamma(st)
        @test γ isa AbstractVector
        @test st.gamma isa AbstractVector
        @test st.gamma_s(st.s) isa AbstractVector
        @test st.gamma_max isa RealU
        @test st.gamma_max_abs isa Real
        @test st.orbit_short isa RealU
        @test st.dsd isa RealU
        @test st.dfs isa RealU
        @test st.dso isa RealU
        @test st.dod isa RealU
    end

    if st isa SinoMoj
        @test st.d_moj(0) isa RealU
        @test st.d_ang isa AbstractVector # angular dependent d for :moj
    end

    if st isa SinoPar
        rs = @inferred rays(st)
        @test rs isa Iterators.ProductIterator
    elseif st isa SinoMoj
        rs = @inferred rays(st)
        @test rs isa Base.Generator{<:Base.Iterators.ProductIterator}
    else
        rs = rays(st) # @NOTinferred because of "fun" closure
        @test rs isa Base.Generator{<:Base.Iterators.ProductIterator}
    end

    @test collect(rs) isa Array{<:Tuple} # @NOTinferred

    x = (1:4) * oneunit(d)
    @inferred st.taufun(x, 2*x)

    ig = ImageGeom( (5,7), (1,1) .* d)
    @test (@inferred footprint_size(st, ig)) isa Float32

    true
end


@testset "orbit-start" begin
    d, orbit = 2mm, 180.0° # stress
    orbit_start = 20 * oneunit(orbit)
    for geo in list
        st = @inferred geo( ; nb = 6, na = 5, d, orbit, orbit_start)
        @test _test_prop(st)
    end
end

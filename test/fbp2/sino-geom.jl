#=
test/fbp2/sino-geom.jl
=#

include("../helper.jl")

using Sinograms: RealU
using Sinograms: SinoGeom, SinoParallel, SinoFan
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
#using Sinograms #: sino_geom, sino_geom_par, sino_geom_fan, sino_geom_moj
using Sinograms: ones, zeros, angles, rays, downsample, oversample # values
using Sinograms: sino_w, sino_s
import Sinograms as SG
# todo

using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred

#=
sg = SinoPar(; d=2mm, orbit=180°)
@inferred sino_w(sg)
@inferred sino_s(sg)
sg = SinoFanFlat(; d=2mm, orbit=180°)
    @inferred SG.sino_geom_gamma(sg)
sg = SinoMoj(; d=2mm, orbit=180°)
@inferred rays(sg)
sg = SinoFanFlat(;d=2.0mm)
@inferred rays(sg)
todo
=#


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
        @inferred SG.sino_geom_taufun(sg, x, 2x)
        @inferred sg.taufun(x, 2*x)
    end
end


# sino_geom.jl

#ig = image_geom(nx=512, fov=500).down(down)
#ig = image_geom(nx=ig.nx, dx=ig.dx, mask=ig.circ())

sg_list = (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)

#=
sg_list = (
    sino_geom(:par ; down=down, d=4),
    sino_geom(:ge1 ; down=down, orbit_start=20, dfs=0), # arc
    sino_geom(:ge1 ; down=down, orbit_start=20, dfs=Inf), # flat
    sino_geom(:moj ; down=down, d=4*sqrt(2)),
    sino_geom(:fan ; down=down, orbit=:short),
)
=#

function _test_prop(sg; d = 2mm, orbit = 180.0°)
    sg.ad[2]
    sg.rfov
    @inferred sg.down(2)
    @inferred sg.over(2)
    sg.dim
    sg.w
    sg.ones
    sg.zeros
    sg.dr
    sg.ds
    sg.r
    sg.s
    sg.ad
    sg.ar
    sg.xds
    sg.yds

    show(isinteractive() ? stdout : devnull, sg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), sg)
    @test (@inferred ones(sg)) isa AbstractMatrix{Float32}
    @test (@inferred zeros(sg)) isa AbstractMatrix{Float32}
    @test (@inferred sino_w(sg)) isa Number
    @test (@inferred sino_s(sg)) isa AbstractVector
    @test (@inferred angles(sg)) isa AbstractVector
    rs = @inferred rays(sg)
    @test rs isa Tuple
    @test (@NOTinferred propertynames(sg)) isa NTuple

    @inferred SG.sino_geom_gamma(sg)

    if sg isa SinoParallel
        @test isnothing(sg.gamma)
    end

    if sg isa SinoFan
        @test sg.gamma isa AbstractVector
        sg.gamma_max
        @test sg.orbit_short isa Real
        sg.dsd
        sg.dfs
        sg.dso
#       sg.dod
    end

    if sg isa SinoMoj
        sg.d_moj(0)
        sg.d_ang # angular dependent d for :moj
    end

    @test sg.shape(vec(sg.ones)) == sg.ones
    x = (1:4) * oneunit(d)
    @inferred sg.taufun(x, 2*x)
    sg.unitv()
    sg.unitv(ib=1, ia=2)
    true
end

@testset "orbit-start" begin
    down = 8
    d, orbit = 2mm, 180.0°
    orbit_start = 20°
    for geo in sg_list
        sg = @inferred geo( ; d, orbit, orbit_start, down)
        @test _test_prop(sg)
    end
end


#=
    sg.plot_grid(plot)

nsg = length(sg_list)
pl = Array{Any}(undef, nsg)
for ii = 1:nsg

# pl[ii] = plot(); sg.plot!(plot! ; ig=ig)

#@inferred todo
@test sino_geom(:ge1 ; units=:cm) isa SinoFanArc

@test_throws String sino_geom(:badhow)
@test_throws String sino_geom(:ge1 ; dfs=-1)
@test_throws String sino_geom(:ge1 ; units=:bad)

pg = sino_geom_plot_grids(plot)
ps = sino_geom_show()

plot(pg..., ps..., layout=(2,4))
prompt()
plot(pl[1:4]...)
=#

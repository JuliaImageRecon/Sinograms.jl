#=
fbp3/ct-geom.jl
CT geometry for 3D tomographic image reconstruction
that describes the sampling characteristics of a CT scan
using the parallel, fan beam arc, or flat fan beam geometry.
2022-05-11, Jason Hu, Jeff Fessler translated from ct_geom.m in MIRT
=#

export dims, downsample, oversample, rays, axes


# Methods common to all types

dims(rg::CtGeom) = (rg.ns, rg.nt, rg.na)

# use ° via `angles()` because mainly for plots
Base.axes(rg::CtGeom) = (_s(rg), _t(rg), angles(rg))

_ws(rg::CtGeom) = (rg.ns-1) // 2 + rg.offset_s
_wt(rg::CtGeom) = (rg.nt-1) // 2 + rg.offset_t

_s(rg::CtGeom) = rg.ds * ((0:rg.ns-1) .- _ws(rg))
_t(rg::CtGeom) = rg.dt * ((0:rg.nt-1) .- _wt(rg))


# down/up sampling

function _downsample(rg::CtGeom, down_s::Int, down_t::Int, down_a::Int)
    ns = 4 * max(rg.ns ÷ 4down_s, 1)
    nt = 2 * max(rg.nt ÷ 2down_t, 1)
    na = max(rg.na ÷ down_t, 1)

    out = (ns, nt, rg.ds * down_s, rg.dt * down_t, rg.offset_s, rg.offset_t,
         na, rg.orbit, rg.orbit_start,
    )
    if rg isa CtFan
         out = (out..., rg.source_offset, rg.dsd, rg.dod)
    end
    return (out..., rg.src)
end


"""
    downsample(rg, down::Int)
    downsample(rg, down::NTuple{3,Int})
Down-sample CT geometry
(for testing with small problems).
"""
downsample(rg::CtGeom, down::Int) = downsample(rg, (down,down,down))

function downsample(rg::G, down::NTuple{3,Real}) where {G <: CtGeom}
    return all(==(1), down) ? rg : G(_downsample(rg, down...)...)
end


function _oversample(rg::CtGeom, over::Int)
    return (
        rg.ns * over, rg.nt * over,
        rg.ds / over, rg.dt / over,
        rg.offset_s * over, rg.offset_t * over,
        rg.na, rg.orbit, rg.orbit_start,
    )
end

"""
    oversample(rg, over::Int)

Over-sample CT geometry in "radial" dimension.
"""
function oversample(rg::G, over::Int) where {G <: CtParallel}
    return (over == 1) ? rg : G(_oversample(rg, over)..., rg.src)
end

function oversample(rg::G, over::Int) where {G <: CtFan}
    return (over == 1) ? rg :
        G(_oversample(rg, over)..., rg.source_offset, rg.dsd, rg.dod, rg.src)
end


"""
    i = rays(rg::CtGeom)
Return parallel-beam coordinates of all rays for this CT geometry.
Return type of `i` is a `ProductIterator` that makes tuples of the form
`(u, v, ϕ, θ)`.
To make projections call
`p = [fun(c...) for c in i]` where `fun` is `radon(...)`.
"""
function rays(rg::CtPar)
    u = _s(rg)
    v = _t(rg)
    ϕ = _ar(rg)
    θ = zero(eltype(ϕ))
    i = Iterators.product(u, v, ϕ, θ)
    return i
end


function rays(rg::CtFan{Td,To}) where {Td,To}
    rg.src isa CtSourceCircle || throw("non-circular not done")
    s = _s(rg)
    t = _t(rg)
    β = _ar(rg)
    i = Iterators.product(s, t, β)
    if rg isa CtFanArc
        fun = stb -> cb_arc_to_par(stb..., _dso(rg), rg.dod)
    else
        fun = stb -> cb_flat_to_par(stb..., _dso(rg), rg.dod)
    end
    return Iterators.map(fun, i)
end



# basic 3d methods


#_cone_angle(rg::CtParallel) = 0
_cone_angle(rg::CtFan) = atan((rg.nt * rg.dt)/2 / rg.dsd)

_zfov(rg::CtParallel) = rg.nt * rg.dt
_zfov(rg::CtFan) = _dso(rg) / rg.dsd * rg.nt * rg.dt


_source_dz_per_view(src::CtSource, na, orbit, zfov::Td) where {Td <: RealU} = zero(Td)
function _source_dz_per_view(src::CtSourceHelix, na, orbit, zfov::Td) where {Td <: RealU}
    na_per_360 = na * (360 / orbit)
    return na == 1 ? zero(Td) : src.pitch * zfov / na_per_360
end

_source_dz_per_view(rg::CtGeom) =
    _source_dz_per_view(rg.src, rg.na, rg.orbit, _zfov(rg))


_source_zs(rg::CtGeom{Td,To,<:CtSourceCircle}) where {Td,To} = fill(zero(Td), rg.na)
#_source_zs(rg::CtGeom{Td,To,<:CtSourceUser}) where {Td,To} = rg.src.source_zs

function _source_zs(rg::CtGeom{Td,To,<:CtSourceHelix}) where {Td,To}
    source_dz = _source_dz_per_view(rg)
    return (0:rg.na-1) * source_dz .+ rg.src.source_z0
end

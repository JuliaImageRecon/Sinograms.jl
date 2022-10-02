#=
fbp3/ct-geom.jl
CT geometry for 3D image reconstruction
that describes the sampling characteristics of a CT scan
using the parallel, fan beam arc, or flat fan beam geometry.
Use sino_geom() instead for the 2D version.

todo: isolate the source trajectory? (e.g., helix, etc.)

2022-05-11, Jason Hu, Jeff Fessler translated from ct_geom.m in MIRT
=#

export CtGeom, CtParallel, CtFan
export CtPar, CtFanArc, CtFanFlat
export ct_geom_ge1
export dims

import Base: ones, zeros

using ImageGeoms: ImageGeom, fovs


"""
    CtGeom{Td,To}

Abstract type for representing 3D CT imaging geometries.
Currently supports only constant pitch source trajectories.

# Common struct elements
* `ns` size of each projection view
* `nt`
* `na` # of projection views
* `orbit` source orbit in degrees
* `orbit_start` starting angle in degrees
* `ds` detector pixel spacing
* `dt`
* `offset_s::Float32` # detector center offset (usually 0)
* `offset_t::Float32`

# Additional elements for `CtFan` types:
* `dsd` distance from source to detector
* `dso` distance from source to origin
* `dod` distance from origin to detector
* `pitch` helix pitch
* `source_z0` initial z position of x-ray source

# Basic methods
* `dims(cg)`
todo

#  Derived values (available by `getproperty`), i.e., `cg.?`

* `.s` s sample locations
* `.t` t sample locations
* `.gamma` [nb] gamma sample values [radians]
* `.gamma_s` gamma values given s values
* `.gamma_max` half of fan angle [radians], if offset_s=0
* `.gamma_max_abs` half of fan angle [radians] - recommended
* `.ws` (ns-1)/2 + st.offset_s
* `.wt` (nt-1)/2 + st.offset_t
* `.ad` [na] source angles in degrees
* `.ar` [na] source angles in radians
* `.dim` dimensions: [st.ns st.nt st.na]
* `.downsample(down)` reduce sampling by integer factor
* `.ones` ones(ns,nt,na, 'single')
* `.zeros` zeros(ns,nt,na, 'single')
* `.rad` [ns] radial distance of each ray from origin
* `.rfov` max radius within FOV
* `.footprint_size(ig)` max footprint width in 's'
* `.cone_angle` (half) cone angle on axis (s=0): +/- angle
* `.zfov` axial FOV
* `.source_zs` [na] z-locations of source for each view

# Methods
* `.down(down)` reduce sampling by integer factor
* `.shape(proj)` reshape `proj`
* `.plot2(ig)` plot the 2D geometry of the image
* `.plot3(ig)` plot the 3D geometry of the cone beam setup
* `.unitv(is,it,ia)` unit vector with single nonzero element
"""
abstract type CtGeom{Td,To} end

abstract type CtParallel{Td,To} <: CtGeom{Td,To} end #line 238 of ct_geom
abstract type CtFan{Td,To} <: CtGeom{Td,To} end


# types


struct CtPar{Td,To} <: CtParallel{Td,To}
    ns::Int
    nt::Int
    na::Int
    orbit::To
    orbit_start::To
    ds::Td
    dt::Td
    offset_s::Float32
    offset_t::Float32
end


struct CtFanArc{Td,To} <: CtFan{Td,To}
    ns::Int
    nt::Int
    na::Int
    orbit::To
    orbit_start::To
    ds::Td
    dt::Td
    offset_s::Float32
    offset_t::Float32
    dsd::Td
#   dso::Td
    dod::Td
#   dfs::Td
    pitch::Float32
#   user_source_zs::Vector{Td}
    source_z0::Td
end


struct CtFanFlat{Td,To} <: CtFan{Td,To}
    ns::Int
    nt::Int
    na::Int
    orbit::To
    orbit_start::To
    ds::Td
    dt::Td
    offset_s::Float32
    offset_t::Float32
    dsd::Td
#   dso::Td
    dod::Td
#   dfs::Td
    pitch::Float32
#   user_source_zs::Vector{Td}
    source_z0::Td
end


# constructors


function CtPar( ;
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
#   orbit::Union{Symbol,Real} = 360,
    orbit::RealU = 360,
    orbit_start::RealU = zero(orbit),
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    down::Int = 1,
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt)
    st = CtPar(ns, nt, na,
            To(orbit), To(orbit_start),
            Td(ds), Td(dt), Float32(offset_s), Float32(offset_t),
    )
    if down > 1
        st = downsample(st, down)
    end
    return st
end


function CtFanArc( ;
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
#   orbit::Union{Symbol,Real} = 360,
    orbit::RealU = 360,
    orbit_start::RealU = zero(orbit),
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    dsd::RealU = 4ns * oneunit(ds),
    dod::RealU = zero(ds),
    down::Int = 1,
    pitch::Real = 0,
    source_z0::RealU = zero(ds),
#   user_source_zs = [],
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt, dsd, dod, source_z0)
    st = CtFanArc(ns, nt, na,
            To(orbit), To(orbit_start),
            Td(ds), Td(dt), Float32(offset_s), Float32(offset_t), Td(dsd), Td(dod),
            Float32(pitch), Td(source_z0), # Td.(user_source_zs,
    )
    if down > 1
        st = downsample(st, down)
    end
    return st
end


function CtFanFlat( ;
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
#   orbit::Union{Symbol,Real} = 360,
    orbit::RealU = 360,
    orbit_start::RealU = zero(orbit),
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    dsd::RealU = 4ns * oneunit(ds),
    dod::RealU = zero(ds),
    down::Int = 1,
    pitch::Real = 0,
    source_z0::RealU = zero(ds),
#   source_zs = (0:st.na-1) * source_dz .+ source_z0
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt, dsd, dod, source_z0)
    st = CtFanFlat(ns, nt, na,
            To(orbit), To(orbit_start),
            Td(ds), Td(dt), Float32(offset_s), Float32(offset_t), Td(dsd), Td(dod),
            Float32(pitch), Td(source_z0), # Td.(user_source_zs,
    )
    if down > 1
        st = downsample(st, down)
    end
    return st
end



"""
    ct_geom_ge1( ; na::Int = 1)
CT geometry for GE Lightspeed system.
These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl.
"""
function ct_geom_ge1( ; na::Int = 1)
    return CtFanArc( ; ns = 888, nt = 64, na = 984, ds = 1.0239, dt = 1.0964,
        offset_s = 1.25, dsd = 949.075, dod = 408.075,
)
end


# methods that are common to all types

dims(st::CtGeom) = (st.ns, st.nt, st.na)
# angles() (in sino-geom.jl)
orbit_short(st::CtGeom) = 180 + 2 * rad2deg(st.gamma_max)
Base.ones(st::CtGeom) = ones(Float32, dims(st))::Array{Float32,3}
Base.zeros(st::CtGeom) = zeros(Float32, dims(st))::Array{Float32,3}

ct_geom_ws(st::CtGeom) = ((st.ns-1)//2 + st.offset_s)::Toffset
ct_geom_wt(st::CtGeom) = ((st.nt-1)//2 + st.offset_t)::Toffset
#ct_geom_s(st::CtGeom) = st.ds * ((0:st.ns-1) .- st.ws) # can't infer
function ct_geom_s(
    st::CtGeom{Td} ;
    T::DataType = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.ds, st.ws, st.ns)
end
function ct_geom_t(
    st::CtGeom{Td} ;
    T::DataType = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.dt, st.wt, st.nt)
end

# basic methods

# just to make the other parts of the code more convenient that use dfs
ct_geom_dfs(st::CtFanArc{Td}) where Td = zero(Td)
ct_geom_dfs(st::CtFanFlat{Td}) where Td = Inf * oneunit(Td)

ct_geom_rfov(st::CtPar) = maximum(abs, st.s)
ct_geom_rfov(st::CtFanArc) = st.dso * sin(ct_geom_gamma_max_abs(st))
ct_geom_rfov(st::CtFanFlat) = st.dso * sin(ct_geom_gamma_max_abs(st)) # todo: correct?

ct_geom_xds(st::CtFanArc) = st.dsd * sin.(st.gamma_s(st.s))
ct_geom_xds(st::CtFanFlat) = st.s
ct_geom_xds(st::CtParallel) = st.s

ct_geom_yds(st::CtFanArc) = @. st.dso - st.dsd * cos(st.gamma_s(st.s))
ct_geom_yds(st::CtFanFlat) = @. fill(-st.dod, st.ns)
ct_geom_yds(st::CtParallel) = 0 * st.s

ct_geom_cone_angle(st::CtParallel) = 0
ct_geom_cone_angle(st::CtFan) = atan((st.nt * st.dt)/2 / st.dsd)

ct_geom_zfov(st::CtParallel) = st.nt * st.dt
ct_geom_zfov(st::CtFan) = st.dso / st.dsd * st.nt * st.dt

ct_geom_gamma_s(st::CtFanArc, ss) = ss / st.dsd
ct_geom_gamma_s(st::CtFanFlat, ss) = atan.(ss / st.dsd)

ct_geom_gamma(st::CtGeom) = ct_geom_gamma_s(st, st.s)
ct_geom_gamma_max(st::CtGeom) = maximum(st.gamma)
ct_geom_gamma_max_abs(st::CtGeom) = maximum(abs, st.gamma)


# type inference help:
_rays_type(Td,To) = Iterators.ProductIterator{
    Tuple{LinRange{Td, Int}, LinRange{Td, Int}, LinRange{To, Int}, To}
}

"""
    i = rays(st::CtGeom)
Return parallel-beam coordinates of all rays for this CT geometry.
Return type of `i` is a `ProductIterator` that makes tuples of the form
`(u, v, ϕ, θ)`.
To make projections call
`p = [fun(c...) for c in i]` where `fun` is `radon(...)`.
"""
function rays(st::CtPar{Td,To})::_rays_type(Td,To) where {Td,To}
    u = st.s
    v = st.t
    ϕ = st.ar / oneunit(eltype(st.ar))
    θ = zero(eltype(ϕ))
    i = Iterators.product(u, v, ϕ, θ)
    return i
end


function rays(st::CtFan{Td,To}) where {Td,To}
    st.pitch == 0 || throw("pitch not done")
    s = st.s
    t = st.t
    β = st.ar / oneunit(eltype(st.ar))
    i = Iterators.product(s, t, β)
    if st isa CtFanArc
        fun = stb -> cb_arc_to_par(stb..., st.dso, st.dod)
    else
        fun = stb -> cb_flat_to_par(stb..., st.dso, st.dod)
    end
    return Iterators.map(fun, i)
end


function ct_geom_source_dz_per_view(st::CtGeom{Td}) where Td
    if st.na == 1 || st.pitch == 0
        return zero(Td)
    end
    na_per_360 = st.na * (360 / st.orbit)
    out = st.pitch * st.zfov / na_per_360
end

function ct_geom_source_zs(st::CtGeom)
#   if !isempty(st.user_source_zs)
#       return st.user_source_zs
#   else
        source_dz = st.source_dz_per_view
        return (0:st.na-1) * source_dz .+ st.source_z0
#   end
end


function ct_geom_footprint_size(st::CtPar, ig::ImageGeom{3})::Float32
    di = sqrt(sum(abs2, ig.deltas[1:2]))
    smax = maximum(abs, st.s)
    return di / st.ds
end

function ct_geom_footprint_size(st::CtFanFlat, ig::ImageGeom{3})::Float32
    di = sqrt(sum(abs2, ig.deltas[1:2]))
    smax = maximum(abs, st.s)
    rfov = maximum(fovs(ig)[1:2]) / 2
    rfov > 0.99 * st.dso && throw("bad dso")
    return di / st.ds * sqrt(st.dsd^2 + smax^2) / (st.dso - rfov)
end

function ct_geom_footprint_size(st::CtFanArc, ig::ImageGeom{3})::Float32
    di = sqrt(sum(abs2, ig.deltas[1:2]))
    smax = maximum(abs, st.s)
    rfov = maximum(fovs(ig)[1:2]) / 2
    rfov > 0.99 * st.dso && throw("bad dso")
    return di / st.ds * st.dsd / (st.dso - rfov)
end


function _downsample(st::CtGeom, down_s::Int, down_t::Int, down_a::Int)
    ns = 4 * max(st.ns ÷ 4down_s, 1)
    nt = 2 * max(st.nt ÷ 2down_t, 1)
    na = max(st.na ÷ down_t, 1)

#=
    user_spec = false
    if !isempty(st.user_source_zs)
        st.user_source_zs = st.user_source_zs[1:down_a:st.na]
        user_spec = true
    end

    if length(st.orbit_start) > 1
        st.orbit_start = st.orbit_start[1:down_a:si.na]
        user_spec = true
    end
=#

#=
    if user_spec
        st.na = length([1:down_a:st.na])
    elseif all(diff(si.source_zs) == 0)
        st.na = max(1, round(si.na / down_a))
    else
        nturn = si.orbit / 360
        na1 = si.na / nturn

        tol = 0.0001
        if abs(round(na1)-na1) < tol
            na2 = round(na1 / down_a)
            tmp = nturn * na2
            if abs(round(tmp)-tmp) < tol
                st.na = round(tmp)
            else
                st.na = floor(tmp)
                st.orbit = 360 / na2 * st.na
            end
        else
            st.na = round(si.na / down_a)
        end
    end
=#

    out = (ns, nt, na, st.orbit, st.orbit_start,
         st.ds * down_s, st.dt * down_t, st.offset_s, st.offset_t)
    if st isa CtFan
         out = (out..., st.dsd, st.dod, st.pitch, st.source_z0)
    end
    return out
end


"""
    downsample(ct, down::Int)
    downsample(ct, down::NTuple{3,Int})
Down-sample CT geometry (for testing with small problems).
"""
downsample(ct::CtGeom, down::Int) = downsample(ct, (down,down,down))

function downsample(ct::C, down::NTuple{3,Real}) where {C <: CtGeom}
    return all(==(1), down) ? ct : C(_downsample(ct, down...)...)::C
end

function ct_geom_unitv(st::CtGeom ;
    is::Int = floor(Int, st.ns/2) + 1,
    it::Int = floor(Int, st.nt/2) + 1,
    ia::Int = 1,
)
    out = zeros(st)
    out[is, it, ia] = 1
    return out
end


ct_geom_fun0 = Dict([
#   (:help, st -> ct_geom_help()),

#   (:dim, st -> (st.ns, st.nt, st.na)),
    (:ws, st -> ct_geom_ws(st)),
    (:wt, st -> ct_geom_wt(st)),
#   (:ones, st -> ones(Float64, st.dim)),
#   (:zeros, st -> zeros(Float64, st.dim)),
    (:unitv, st -> ((; kwarg...) -> ct_geom_unitv(st ; kwarg...))),

    (:s, st -> ct_geom_s(st)),
    (:t, st -> ct_geom_t(st)),
    (:dso, st -> st.dsd - st.dod),
    (:dfs, st -> ct_geom_dfs(st)),

    (:gamma, st -> ct_geom_gamma(st)),
    (:gamma_s, st -> (ss -> ct_geom_gamma_s(st, ss))),
    (:gamma_max, st -> ct_geom_gamma_max(st)),
    (:gamma_max_abs, st -> ct_geom_gamma_max_abs(st)),
    (:xds, st -> ct_geom_xds(st)),
    (:yds, st -> ct_geom_yds(st)),

    (:cone_angle, st -> ct_geom_cone_angle(st)),
    (:zfov, st -> ct_geom_zfov(st)),
    (:source_dz_per_view, st -> ct_geom_source_dz_per_view(st)),
    (:source_zs, st -> ct_geom_source_zs(st)),
    (:orbit_short, st -> orbit_short(st)),
#   (:downsample, st -> (down::Int -> ct_geom_downsample(st, down))),

    (:ad, st -> angles(st)),
    (:ar, st -> to_radians(st.ad)),
#   (:rad, st -> ct_geom_rad(st)),
    (:rfov, st -> ct_geom_rfov(st)),
    (:shape, st -> (proj -> reshaper(proj, dims(st)))),
    (:footprint_size, st -> (ig -> ct_geom_footprint_size(st, ig))),

#   (:plot_grid, st -> ((plot::Function) -> ct_geom_plot_grid(st, plot))),
#   (:plot2, st -> ((ig::ImageGeom) -> ct_geom_plot2(st, ig))),
#   (:plot3, st -> ((ig::ImageGeom) -> ct_geom_plot3(st, ig)))
])

Base.getproperty(st::CtGeom, s::Symbol) =
    haskey(ct_geom_fun0, s) ? ct_geom_fun0[s](st) :
    getfield(st, s)

Base.propertynames(st::CtGeom) =
    (fieldnames(typeof(st))..., keys(ct_geom_fun0)...)

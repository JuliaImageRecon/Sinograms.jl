#=
geom/type3.jl
CT geometry for 3D tomographic image reconstruction
that describes the sampling characteristics of a CT scan
using the parallel, fan beam arc, or flat fan beam geometry.

2022-05-11, Jason Hu, Jeff Fessler translated from ct_geom.m in MIRT
=#

export CtGeom, CtParallel, CtFan
export CtPar, CtFanArc, CtFanFlat


"""
    CtGeom{Td,To,Ts}

Abstract type for representing ray geometries
for 3D CT imaging.

The projection view coordinates are `(s,t)`
where `s` denotes the transaxial sampling
and `t` denotes the axial direction (along `z`).

# Common fields

* `ns` size of each projection view
* `nt`
* `ds` detector pixel spacing
* `dt`
* `offset_s` unitless detector center offset (usually 0 or 0.25)
* `offset_t` unitless, usually 0
* `na` # of projection views, aka nϕ or nβ
* `orbit` source orbit in degrees (or Unitful)
* `orbit_start` starting angle in degrees (or Unitful)
  `orbit` and `orbit_start` must both be unitless (degrees) or have same units.
* `src::CtSource` describes the X-ray CT source trajectory.
  Primary support for `CtSourceCircle()`.

# Additional fields for `CtFan` types:

* `dsd` distance from source to detector
* `dod` distance from origin to detector
* `source_offset` usually 0

# Units:

* `ds`, `dt`, `source_offset`, `dsd`, `dod`
  must all be unitless or have the same units.

# Basic methods

* `angles` (na) in degrees
* `dims (ns, nt, na)`
* `ones = ones(Float32, ns,nt,na)`
* `zeros = zeros(Float32, ns,nt,na)`
* `rays` iterator of `(u,v,ϕ,θ)` samples
* `downsample(st, down)` reduce sampling by integer factor
* `oversample(st, over)`
* `ct_geom_plot3` plot system geometry

# Non-exported helper functions for developers:

* `_s (ns) s` sample locations
* `_t (nt) t` sample locations
* `_ws = (ns-1)/2 + offset_s` "middle" sample position
* `_wt = (nt-1)/2 + offset_t`
* `_ar (na)` source angles [radians]
* `_rfov` max radius within FOV
* `_zfov` axial FOV
* `_xds (nb)` center of detector elements (beta=0)
* `_yds (nb)` ""
* `_tau(rg, x, y)` projected s/ds for each `(x,y)` pair `(length(x), na)`
* `_shape(rg, proj [,:])` reshape `proj` into array `(ns,nt,na[,:])`
* `_unitv(rg [, (is,it,ia)])`
  unit 'vector' with single nonzero element

For fan beam:

* `_dso = dsd - dod` distance from source to origin (Inf for parallel beam)
* `_dfs` distance from source to detector focal spot
        (0 for 3rd gen CT, `Inf` for flat detectors)
* `_gamma(rg [,s]) (ns)` gamma sample values `radians`, optionally given `s` values
* `_gamma_max = max(|γ|)` half of fan angle `radians`, if `offset_s` == 0
* `_cone_angle` (half) cone angle on axis (s=0): +/- angle

# Notes
Use `sino_geom()` instead for 2D geometries.
"""
abstract type CtGeom{Td <: RealU, To <: RealU, Ts} <: RayGeom{Td,To} end

abstract type CtParallel{Td <: RealU, To <: RealU, Ts} <: CtGeom{Td,To,Ts} end
abstract type CtFan{Td <: RealU, To <: RealU, Ts} <: CtGeom{Td,To,Ts} end


# types


"""
    CtPar{Td,To,Ts}
3D parallel-beam projection geometry
"""
struct CtPar{Td <: RealU, To <: RealU, Ts} <: CtParallel{Td,To,Ts}
    # detector:
    ns::Int
    nt::Int
    ds::Td
    dt::Td
    offset_s::Toffset
    offset_t::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
    src::Ts
end


#=
"""
    CtMoj{Td,To,Ts}
3D Mojette
where `d` means `dx` (square pixel size)
"""
struct CtMoj{Td <: RealU, To <: RealU, Ts} <: CtParallel{Td,To,Ts}
    # detector:
    ns::Int
    nt::Int
    d::Td
    dt::Td
    offset_s::Toffset
    offset_t::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
    src::Ts
end
=#


"""
    CtFanArc{Td,To,Ts}
3D CBCT geometry for arc detector
"""
struct CtFanArc{Td <: RealU, To <: RealU, Ts} <: CtFan{Td,To,Ts}
    # detector:
    ns::Int
    nt::Int
    ds::Td
    dt::Td
    offset_s::Toffset
    offset_t::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
    source_offset::Td
    dsd::Td
    dod::Td
    src::Ts
end


"""
    CtFanFlat{Td,To,Ts}
3D CTCT geometry for flat detector
"""
struct CtFanFlat{Td <: RealU, To <: RealU, Ts} <: CtFan{Td,To,Ts}
    # detector:
    ns::Int
    nt::Int
    ds::Td
    dt::Td
    offset_s::Toffset
    offset_t::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
    source_offset::Td
    dsd::Td
    dod::Td
    src::Ts
end


# constructors


"""
    CtPar( ; ns nt ds dt offset_s offset_t na orbit orbit_start)

Constructor with named keywords.
See `?CtGeom` for documentation.

```jldoctest
julia> CtPar()
CtPar{Float32, Float32, CtSourceCircle} :
 ns::Int64 128
 nt::Int64 64
 ds::Float32 1.0
 dt::Float32 1.0
 offset_s::Float32 0.0
 offset_t::Float32 0.0
 na::Int64 60
 orbit::Float32 180.0
 orbit_start::Float32 0.0
 src::CtSourceCircle CtSourceCircle()
```
"""
function CtPar( ;
    ns::Int = 128,
    nt::Int = 64,
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    na::Int = 60,
    orbit::RealU = 180,
    orbit_start::RealU = zero(eltype(orbit)),
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt)
    return CtPar(ns, nt, Td(ds), Td(dt), Toffset(offset_s), Toffset(offset_t),
        na, To(orbit), To(orbit_start),
        CtSourceCircle(),
    )
end


#=
"""
    CtMoj( ; ns nt ds dt offset_s offset_t na orbit orbit_start)

Constructor with named keywords.
See `?CtGeom` for documentation.

```jldoctest
julia> CtMoj()
CtMoj{Float32, Float32, CtSourceCircle} :
 ns::Int64 128
 nt::Int64 64
 ds::Float32 1.0
 dt::Float32 1.0
 offset_s::Float32 0.0
 offset_t::Float32 0.0
 na::Int64 60
 orbit::Float32 180.0
 orbit_start::Float32 0.0
 src::CtSourceCircle CtSourceCircle()
```
"""
function CtMoj( ;
    ns::Int = 128,
    nt::Int = 64,
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    na::Int = 60,
    orbit::RealU = 180,
    orbit_start::RealU = zero(eltype(orbit)),
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt)
    return CtMoj(ns, nt, Td(ds), Td(dt), Toffset(offset_s), Toffset(offset_t),
        na, To(orbit), To(orbit_start),
        CtSourceCircle(),
    )
end
=#


"""
    CtFanArc( ; ns nt ds dt offset_s offset_t
        na orbit orbit_start
        dsd = 4ns * ds, dod = ns * ds)
    CtFanArc(:short ; ...)

Constructor with named keywords.
See `?CtGeom` for documentation.

* Use `:short` argument to specify a short scan,
  in which case `na` will be scaled down proportionally as well.

```jldoctest
julia> CtFanArc()
CtFanArc{Float32, Float32, CtSourceCircle} :
 ns::Int64 128
 nt::Int64 64
 ds::Float32 1.0
 dt::Float32 1.0
 offset_s::Float32 0.0
 offset_t::Float32 0.0
 na::Int64 64
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float32 0.0
 dsd::Float32 512.0
 dod::Float32 128.0
 src::CtSourceCircle CtSourceCircle()
```
"""
function CtFanArc( ;
    ns::Int = 128,
    nt::Int = 64,
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    na::Int = 64,
    orbit::RealU = 360,
    orbit_start::RealU = zero(eltype(orbit)),
    source_offset::RealU = zero(eltype(ds)),
    dsd::RealU = 4 * ns * ds,
    dod::RealU = ns * ds,
    src::CtSource = CtSourceCircle(),
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt, source_offset, dsd, dod)
    return CtFanArc(ns, nt, Td(ds), Td(dt), Toffset(offset_s), Toffset(offset_t),
        na, To(orbit), To(orbit_start),
        Td(source_offset), Td(dsd), Td(dod),
        src,
    )
end


# :short case
function CtFanArc(orbit::Symbol ;
    na::Int = 128, orbit_start::To = 0f0, kwargs...,
) where {To <: RealU}
    orbit == :short || error("bad orbit $orbit")
    tmp = CtFanArc( ; orbit=To(360), kwargs...)
    os = _orbit_short(tmp)
    na = ceil(Int, na * os / To(360))
    return CtFanArc( ; na, orbit = To(os), kwargs...)
end


"""
    CtFanFlat( ; ns nt ds dt offset_s offset_t
        na orbit orbit_start
        dsd = 4ns * ds, dod = ns * ds)
    CtFanFlat(:short ; ...)

Constructor with named keywords.
See `?CtGeom` for documentation.

* Use `:short` argument to specify a short scan,
  in which case `na` will be scaled down proportionally as well.

```jldoctest
julia> CtFanFlat()
CtFanFlat{Float32, Float32, CtSourceCircle} :
 ns::Int64 128
 nt::Int64 64
 ds::Float32 1.0
 dt::Float32 1.0
 offset_s::Float32 0.0
 offset_t::Float32 0.0
 na::Int64 64
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float32 0.0
 dsd::Float32 512.0
 dod::Float32 128.0
 src::CtSourceCircle CtSourceCircle()
```
"""
function CtFanFlat( ;
    ns::Int = 128,
    nt::Int = 64,
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    na::Int = 64,
    orbit::RealU = 360,
    orbit_start::RealU = zero(eltype(orbit)),
    source_offset::RealU = zero(eltype(ds)),
    dsd::RealU = 4ns * ds,
    dod::RealU = ns * ds,
    src::CtSource = CtSourceCircle(),
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt, source_offset, dsd, dod)
    return CtFanFlat(ns, nt, Td(ds), Td(dt), Toffset(offset_s), Toffset(offset_t),
        na, To(orbit), To(orbit_start),
        Td(source_offset), Td(dsd), Td(dod),
        src,
    )
end


# :short case
function CtFanFlat(orbit::Symbol ;
    na::Int = 128, orbit_start::To = 0f0, kwargs...,
) where {To <: RealU}
    orbit == :short || error("bad orbit $orbit")
    tmp = CtFanFlat( ; orbit=To(360), kwargs...)
    os = _orbit_short(tmp)
    na = ceil(Int, na * os / To(360))
    return CtFanFlat( ; na, orbit = os, kwargs...)
end


"""
    CtFanArc(::Val{:ge1} ; kwargs...)
GE Lightspeed system CT geometry.

# option
* `unit::RealU = 1` or use `1mm`
* (see `CtFanArc`)

# out
* `CtFanArc`

These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl.

```jldoctest
julia> CtFanArc(Val(:ge1))
CtFanArc{Float64, Float32, CtSourceCircle} :
 ns::Int64 888
 nt::Int64 64
 ds::Float64 1.0239
 dt::Float64 1.0964
 offset_s::Float32 1.25
 offset_t::Float32 0.0
 na::Int64 984
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float64 0.0
 dsd::Float64 949.075
 dod::Float64 408.075
 src::CtSourceCircle CtSourceCircle()
```
"""
function CtFanArc(::Val{:ge1} ;
    unit::RealU = 1,
    ns::Int = 888,
    nt::Int = 64,
    ds::RealU = 1.0239 * unit,
    dt::RealU = 1.0964 * unit,
    offset_s::Real = 1.25,
    offset_t::Real = 0,
    na::Int = 984,
    orbit::Union{Symbol,Real} = 360,
    dsd::RealU = 949.075 * unit,
    dod::RealU = 408.075 * unit,
    kwargs...,
)

    if orbit === :short
        na = 642 # trick: reduce na for short scans
        orbit = na / 984 * 360
    end

    return CtFanArc( ; ns, nt, ds, dt, offset_s, offset_t,
        na, orbit, dsd, dod,
        kwargs...,
    )
end

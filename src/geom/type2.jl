#=
geom/type2.jl
sinogram geometry definitions for 2D tomographic image reconstruction
2019-07-01, Jeff Fessler, University of Michigan
2022-01-22, copied from MIRT.jl and updated to support Unitful values
=#

export SinoGeom, SinoParallel, SinoFan
export SinoPar, SinoMoj, SinoFanArc, SinoFanFlat


"""
    SinoGeom{Td,To}

Abstract type for representing ray geometries
for 2D sinograms.
This describes the sampling characteristics of a given sinogram
for a 2D parallel or fan-beam system.

# Common fields

* `nb` # of "radial" samples, aka `nr` or `ns`
* `d` aka `dr` or `ds`, "radial" sample spacing
* `offset` unitless sample offset (usually 0 or 0.25)
* `na` # of angular samples, aka `nϕ` or `nβ`
   default: `2 * floor(Int, nb * π/2 / 2)`
* `orbit` source orbit in degrees (or Unitful)
* `orbit_start` starting angle in degrees (or Unitful)
  `orbit` and `orbit_start` must both be unitless (degrees) or have same units.

# Additional fields for `SinoFan` types

* `dsd` distance from source to detector
* `dod` distance from origin to detector
* `source_offset` usually 0

# Units:

* `d`, `source_offset`, `dsd`, `dod`
  must all have the same units.

# Basic methods

* `angles` (na) in degrees
* `dims (nb, na)`
* `ones = ones(Float32, nb, na)`
* `zeros = zeros(Float32, nb, na)`
* `rays` iterator of `(r, ϕ)` parallel-beam coordinate tuples of size `(nb, na)`
* `downsample(st, down)` reduce sampling by integer factor
* `oversample(st, over)`
* `sino_geom_plot!` plot system geometry

# Non-exported helper functions for developers:

* `_ds|dr` radial sample spacing (`NaN` for `:moj`)
* `_s (nb) s` sample locations
* `_w = (nb-1)/2 + offset` "middle" sample position
* `_ar (na)` source angles [radians]
* `_rfov` radial FOV
* `_xds (nb)` center of detector elements (beta=0)
* `_yds (nb)` ""
* `_tau(rg, x, y)` projected s/ds for each `(x,y)` pair `(length(x), na)`
* `_shape(rg, sino [,:])` reshape `sino` into array `(nb,na[,:])`
* `_unitv(rg [, (ib,ia)])`
  unit 'vector' with single nonzero element

For mojette:

* `_d_ang (na)` angle-dependent radial spacing

For fan beam:

* `_dso = dsd - dod` distance from source to origin (Inf for parallel beam)
* `_dfs` distance from source to detector focal spot
        (0 for 3rd gen CT, `Inf` for flat detectors)
* `_gamma(rg [,s]) (nb)` gamma sample values `radians`, optionally given `s` values
* `_gamma_max = max(|γ|)` half of fan angle `radians`, if `offset_s` == 0

# Notes
* Use `ct_geom()` instead for 3D axial or helical cone-beam CT.
"""
abstract type SinoGeom{Td <: RealU, To <: RealU} <: RayGeom{Td,To} end

abstract type SinoParallel{Td <: RealU, To <: RealU} <: SinoGeom{Td,To} end
abstract type SinoFan{Td <: RealU, To <: RealU} <: SinoGeom{Td,To} end


# types


"""
    SinoPar{Td,To}
2D parallel-beam sinogram geometry
"""
struct SinoPar{Td <: RealU, To <: RealU} <: SinoParallel{Td,To}
    # detector:
    nb::Int
    d::Td
    offset::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
end


"""
    SinoMoj{Td,To}
2D Mojette sinogram geometry
where `d` means `dx` (square pixel size)
"""
struct SinoMoj{Td <: RealU, To <: RealU} <: SinoParallel{Td,To}
    # detector:
    nb::Int
    d::Td
    offset::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
end


"""
    SinoFanArc{Td,To}
2D fan-beam sinogram geometry for arc detector
"""
struct SinoFanArc{Td <: RealU, To <: RealU} <: SinoFan{Td,To}
    # detector:
    nb::Int
    d::Td
    offset::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
    source_offset::Td
    dsd::Td
    dod::Td
end


"""
    SinoFanFlat{Td,To}
2D fan-beam sinogram geometry for flat detector
"""
struct SinoFanFlat{Td <: RealU, To <: RealU} <: SinoFan{Td,To}
    # detector:
    nb::Int
    d::Td
    offset::Toffset
    # source:
    na::Int
    orbit::To
    orbit_start::To
    source_offset::Td
    dsd::Td
    dod::Td
end


# constructors


"""
    SinoPar( ; nb d offset na orbit orbit_start)

Constructor with named keywords.
See `?SinoGeom` for documentation.

```jldoctest
julia> SinoPar()
SinoPar{Float32, Float32} :
 nb::Int64 128
 d::Float32 1.0
 offset::Float32 0.0
 na::Int64 200
 orbit::Float32 180.0
 orbit_start::Float32 0.0
```
"""
function SinoPar( ;
    nb::Int = 128,
    d::RealU = 1f0,
    offset::Real = 0,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    orbit::RealU = 180,
    orbit_start::RealU = zero(typeof(orbit)),
)

    To = _promoter(orbit, orbit_start)
    return SinoPar(nb, d, Toffset(offset),
        na, To(orbit), To(orbit_start),
    )
end


"""
    SinoMoj( ; nb d offset na orbit orbit_start)

Constructor with named keywords.
See `?SinoGeom` for documentation.

```jldoctest
julia> SinoMoj()
SinoMoj{Float32, Float32} :
 nb::Int64 128
 d::Float32 1.0
 offset::Float32 0.0
 na::Int64 200
 orbit::Float32 180.0
 orbit_start::Float32 0.0
```
"""
function SinoMoj( ;
    nb::Int = 128,
    d::RealU = 1f0,
    offset::Real = 0,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    orbit::RealU = 180,
    orbit_start::RealU = zero(typeof(orbit)),
)

    To = _promoter(orbit, orbit_start)
    return SinoMoj(nb, d, Toffset(offset),
        na, To(orbit), To(orbit_start),
    )
end


"""
    SinoFanArc( ; nb d offset na orbit orbit_start
        source_offset, dsd = 4 * nb * d, dod = nb * d)
    SinoFanArc(:short ; ...)

Constructor with named keywords.
See `?SinoGeom` for documentation.

* `d`, `source_offset`, `dsd`, `dod`
  must all have the same units.

* Use `:short` argument to specify a short scan,
  in which case `na` will be scaled down proportionally as well.

```jldoctest
julia> SinoFanArc()
SinoFanArc{Float32, Float32} :
 nb::Int64 128
 d::Float32 1.0
 offset::Float32 0.0
 na::Int64 200
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float32 0.0
 dsd::Float32 512.0
 dod::Float32 128.0
```
"""
function SinoFanArc( ;
    nb::Int = 128,
    d::RealU = 1,
    offset::Real = 0,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    orbit::RealU = 360,
    orbit_start::RealU = zero(typeof(orbit)),
    source_offset::RealU = zero(typeof(d)),
    dsd::RealU = 4 * nb * d,
    dod::RealU = nb * d,
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(d, source_offset, dsd, dod)
    return SinoFanArc(nb, Td(d), Toffset(offset),
        na, To(orbit), To(orbit_start),
        Td(source_offset), Td(dsd), Td(dod),
    )
end


# :short case
function SinoFanArc(orbit::Symbol ;
    na::Int = 128, orbit_start::To = 0f0, kwargs...,
) where {To <: RealU}
    orbit == :short || error("bad orbit $orbit")
    tmp = SinoFanArc( ; orbit=To(360), kwargs...)
    na = ceil(Int, na * _orbit_short(tmp) / To(360))
    return SinoFanArc( ; na, orbit = To(_orbit_short(tmp)), kwargs...)
end


"""
    SinoFanFlat( ; nb d offset na orbit orbit_start
        source_offset, dsd= 4 * nb * d, dod = nb * d)
    SinoFanFlat(:short ; ...)

Constructor with named keywords.
See `?SinoGeom` for documentation.

* Use `:short` argument to specify a short scan,
  in which case `na` will be scaled down proportionally as well.

```jldoctest
julia> SinoFanFlat()
SinoFanFlat{Float32, Float32} :
 nb::Int64 128
 d::Float32 1.0
 offset::Float32 0.0
 na::Int64 200
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float32 0.0
 dsd::Float32 512.0
 dod::Float32 128.0
```
"""
function SinoFanFlat( ;
    nb::Int = 128,
    d::RealU = 1,
    offset::Real = 0,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    orbit::RealU = 360,
    orbit_start::RealU = zero(typeof(orbit)),
    source_offset::RealU = zero(typeof(d)),
    dsd::RealU = 4 * nb * d,
    dod::RealU = nb * d,
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(d, source_offset, dsd, dod)
    return SinoFanFlat(nb, Td(d), Toffset(offset),
        na, To(orbit), To(orbit_start),
        Td(source_offset), Td(dsd), Td(dod),
    )
end


# :short case
function SinoFanFlat(orbit::Symbol ;
    na::Int = 128, orbit_start::To = 0f0, kwargs...,
) where {To <: RealU}
    orbit == :short || error("bad orbit $orbit")
    tmp = SinoFanFlat( ; orbit=To(360), kwargs...)
    na = ceil(Int, na * _orbit_short(tmp) / To(360))
    return SinoFanFlat( ; na, orbit = To(_orbit_short(tmp)), kwargs...)
end


"""
    SinoFanArc(Val(:ge1) ; kwargs...)
GE Lightspeed system CT geometry.

# option
* `unit::RealU = 1` or use `1mm`
* `nb::Int = 888` # of detector channels
* `d::RealU = 1.0239` channel spacing
* `offset::Real = 1.25` for "quarter-detector" offset
* `na::Int = 984` # of angles
* `orbit::Union{Symbol,Real} = 360` use `:short` for short scan
* `dsd::RealU = 949.075`
* `dod::RealU = 408.075`

# out
* `SinoFanArc`

These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl.

```jldoctest
julia> SinoFanArc(Val(:ge1))
SinoFanArc{Float32, Float32} :
 nb::Int64 888
 d::Float32 1.0239
 offset::Float32 1.25
 na::Int64 984
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float32 0.0
 dsd::Float32 949.075
 dod::Float32 408.075
```
"""
function SinoFanArc(::Val{:ge1} ;
    unit::RealU = 1,
    nb::Int = 888,
    d::RealU = 1.0239f0 * unit,
    offset::Real = 1.25f0,
    na::Int = 984,
    orbit::Union{Symbol,Real} = 360,
    dsd::RealU = 949.075f0 * unit,
    dod::RealU = 408.075f0 * unit,
    kwargs...,
)

    if orbit === :short
        na = 642 # trick: reduce na for short scans
        orbit = Float32(na / 984 * 360)
    end

    return SinoFanArc( ; nb, d, offset,
        na, orbit, dsd, dod,
        kwargs...,
    )
end

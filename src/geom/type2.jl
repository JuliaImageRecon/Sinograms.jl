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

* `nb` # of "radial" samples, aka nr or ns
* `d` aka dr or ds, "radial" sample spacing
* `offset` unitless sample offset (usually 0 or 0.25)
* `na` # of angular samples, aka nϕ or nβ
   default: `2 * floor(Int, nb * π/2 / 2)`
* `orbit` source orbit in degrees (or Unitful)
* `orbit_start` starting angle in degrees (or Unitful)
  `orbit` and `orbit_start` must both be unitless (degrees) or have same units.

# Additional fields for `SinoFan` types

* `dsd` distance from source to detector
* `dod` distance from origin to detector

# Units:

* `d`, `source_offset`, `dsd`, `dod`
  must all have the same units.

# Derived values (available by `getproperty`), i.e., `st.?`

* `.dim` dimensions: `(nb,na)`
* `.ds|dr` radial sample spacing (`NaN` for `:moj`)
* `.s (nb) s` sample locations
* `.w = (nb-1)/2 + offset` "middle" sample position
* `.ad (na)` source angles [degrees]
* `.ar (na)` source angles [radians]
* `.ones`      `ones(Float32, nb,na)`
* `.zeros`     `zeros(Float32, nb,na)`
* `.rfov` radial FOV
* `.xds (nb)` center of detector elements (beta=0)
* `.yds (nb)` ""
* `.plot_grid(scatter)` plot `sg.grid` using `Plots.scatter`

For mojette:

* `.d_ang (na)` angle-dependent radial spacing

For fan beam:

* `source_offset:` same units as d, etc., e.g., [mm] (use with caution!)
* `dso = dsd - dod` distance from source to origin (Inf for parallel beam)
* `dfs` distance from source to detector focal spot
        (0 for 3rd gen CT, `Inf` for flat detectors)
* `.gamma (nb)` gamma sample values [radians]
* `.gamma_max` half of fan angle [radians], if offset=0

# Basic methods

* `dims` (nb, na)
* `sino_s` [nb] s sample locations
* `sino_w` `(nb-1)/2 + offset` ('middle' sample position)
* `rays` `(rgrid, ϕgrid)` [nb na] parallel-beam coordinates
* `ones` `ones(Float32, nb, na)`
* `zeros` `zeros(Float32, nb, na)`
* `downsample`
* `oversample`
* `angles` (na) in degrees

# Methods

* `.shape(sino)` reshape sinograms into array [nb na :]
* `.unitv(;ib,ia)` unit 'vector' with single nonzero element
* `.taufun(x,y)` projected s/ds for each (x,y) pair [numel(x) na]

# Notes
* Use `ct_geom()` instead for 3D axial or helical cone-beam CT.
"""
abstract type SinoGeom{Td,To} <: RayGeom{Td,To} end

abstract type SinoParallel{Td,To} <: SinoGeom{Td,To} end
abstract type SinoFan{Td,To} <: SinoGeom{Td,To} end


# types


"""
    SinoPar
2D parallel-beam sinogram geometry
"""
struct SinoPar{Td,To} <: SinoParallel{Td,To}
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
    SinoMoj
2D Mojette sinogram geometry
where `d` means `dx` (square pixel size)
"""
struct SinoMoj{Td,To} <: SinoParallel{Td,To}
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
   SinoFanArc
2D fan-beam sinogram geometry for arc detector
"""
struct SinoFanArc{Td,To} <: SinoFan{Td,To}
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
    SinoFanFlat
2D fan-beam sinogram geometry for flat detector
"""
struct SinoFanFlat{Td,To} <: SinoFan{Td,To}
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
    orbit_start::RealU = zero(eltype(orbit)),
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
    orbit_start::RealU = zero(eltype(orbit)),
)

    To = _promoter(orbit, orbit_start)
    return SinoMoj(nb, d, Toffset(offset),
        na, To(orbit), To(orbit_start),
    )
end


"""
    SinoFanArc( ; nb d offset na orbit orbit_start
        source_offset, dsd = 4 * nb * d, dod = nb * d)

Constructor with named keywords.
See `?SinoGeom` for documentation.

* `d`, `source_offset`, `dsd`, `dod`
  must all have the same units.

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
    orbit_start::RealU = zero(eltype(orbit)),
    source_offset::RealU = zero(eltype(d)),
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


"""
    SinoFanFlat( ; nb d offset na orbit orbit_start
        source_offset, dsd= 4 * nb * d, dod = nb * d)

Constructor with named keywords.
See `?SinoGeom` for documentation.

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
    orbit_start::RealU = zero(eltype(orbit)),
    source_offset::RealU = zero(eltype(d)),
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


"""
    SinoFan(Val(:ge1) ; kwargs...)
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
SinoFanArc{Float64, Float32} :
 nb::Int64 888
 d::Float64 1.0239
 offset::Float32 1.25
 na::Int64 984
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 source_offset::Float64 0.0
 dsd::Float64 949.075
 dod::Float64 408.075
```
"""
function SinoFanArc(::Val{:ge1} ;
    unit::RealU = 1,
    nb::Int = 888,
    d::RealU = 1.0239 * unit,
    offset::Real = 1.25,
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

    return SinoFanArc( ; nb, d, offset,
        na, orbit, dsd, dod,
        kwargs...,
    )
end
#=
geom/sino-geom.jl
sinogram geometry definitions for 2D tomographic image reconstruction
2019-07-01, Jeff Fessler, University of Michigan
2022-01-22, copied from MIRT.jl and updated to support Unitful values
=#

export dims, downsample, oversample, rays, axes


# Methods common to all types

dims(rg::SinoGeom) = (rg.nb, rg.na)

# use ° via `angles()` because mainly for plots
Base.axes(rg::SinoGeom) = (_s(rg), angles(rg))

_w(rg::SinoGeom) = (rg.nb-1)/2 + rg.offset

_s(rg::SinoGeom) = rg.d * ((0:rg.nb-1) .- _w(rg))
#_r = _s

# This helps unify SinoGeom and CtGeom cases:
_ds(rg::SinoGeom) = rg.d
#_ds(::SinoMoj) = error("undefined")
_ds(rg::CtGeom) = rg.ds
#_ds(::CtMoj) = error("undefined")


# down/up sampling

function _downsample(rg::SinoGeom, down::Int)
    nb = 2 * max(rg.nb ÷ 2down, 1) # keep it even
    na = max(rg.na ÷ down, 1)
    return (nb, rg.d * down, rg.offset, na, rg.orbit, rg.orbit_start)
end


"""
    downsample(rg, down)

Down-sample sinogram geometry
(for testing with small problems).
"""
function downsample(rg::G, down::Int) where {G <: SinoParallel}
    return (down == 1) ? rg : G(_downsample(rg, down)...)
end

function downsample(rg::G, down::Int) where {G <: SinoFan}
    return (down == 1) ? rg :
        G(_downsample(rg, down)..., rg.source_offset, rg.dsd, rg.dod)
end


function _oversample(rg::SinoGeom, over::Int)
    return (
        rg.nb * over, rg.d / over, rg.offset * over,
        rg.na, rg.orbit, rg.orbit_start,
    )
end

"""
    oversample(rg, over::Int)

Over-sample sinogram geometry in "radial" dimension.
For Mojette sampling, it means that `d = dx/over`.
"""
function oversample(rg::G, over::Int) where {G <: SinoParallel}
    return (over == 1) ? rg : G(_oversample(rg, over)...)
end

function oversample(rg::G, over::Int) where {G <: SinoFan}
    return (over == 1) ? rg :
        G(_oversample(rg, over)..., rg.source_offset, rg.dsd, rg.dod)
end


#=
# gamma for general finite dfs (unsupported, would need new SinoFan subtype with dfs)
function _gamma(rg::SinoFan)
    dis_foc_det = _dfs(rg) + rg.dsd
    α = _s(rg) / dis_foc_det # equivalent to s/dsd when dfs=0
    return atan.(dis_foc_det * sin.(α), dis_foc_det * cos.(α) .- _dfs(rg))
end
=#


"""
    i = rays(rg::SinoGeom)

Radial `r` and angular `ϕ` coordinates (in radians)
of all sinogram elements
for the given geometry.
Return type of `i` is a `ProductIterator` that makes tuples of the form
`(r, ϕ)`.
To make projections call
`p = [fun(c...) for c in i]` where `fun` is `radon(...)`.
"""
function rays(rg::SinoPar)
    s = _s(rg)
    ϕ = _ar(rg)
    i = Iterators.product(s, ϕ)
    return i
end


_moj_to_par(s, ϕ) = (s * maximum(abs, sincos(ϕ)), ϕ)
_moj_to_par(sϕ::Tuple) = _moj_to_par(sϕ...)

function rays(rg::SinoMoj)
    s = _s(rg)
    ϕ = _ar(rg)
    i = Iterators.product(s, ϕ)
    return Iterators.map(_moj_to_par, i)
end


function _fan_to_par(rg::SinoFan, s, β)
    γ = _gamma(rg, s)
    r = _dso(rg) * sin(γ)
    ϕ = γ + β
    return (r, ϕ)
end
_fan_to_par(rg::SinoFan, sβ::Tuple) = _fan_to_par(rg, sβ...)

function rays(rg::SinoFan)
    s = _s(rg)
    β = _ar(rg)
    i = Iterators.product(s, β)
    fun = sβ -> _fan_to_par(rg, sβ)
    return Iterators.map(fun, i)
end

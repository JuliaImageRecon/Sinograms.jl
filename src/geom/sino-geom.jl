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
_r = _s

#_dr(rg::SinoGeom) = rg.d
#_ds(rg::SinoGeom) = rg.d
#_dr(::SinoMoj) = NaN
#_ds(::SinoMoj) = NaN


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
    fun = sβ -> _fan_to_par(rg, sβ) # closure prevents type inference? todo
    return Iterators.map(fun, i)
end



# τ (unitless)
function _tau(ϕ::RealU, x::RealU, y::RealU, dr::RealU)
    sϕ, cϕ = sincos(ϕ)
    return (x * cϕ + y * sϕ) / dr
end

function _tau(
    rg::Union{SinoPar,CtPar},
    x::AbstractVector,
    y::AbstractVector,
)
    return _tau.(_ar(rg)', x, y, rg isa SinoPar ? rg.d : rg.ds) # outer-product
end

# this one may not be useful but it helps unify tests
function _tau(rg::SinoMoj, x::AbstractVector, y::AbstractVector)
    d_moj = _d_moj(rg)
    return _tau.(_ar(rg)', x, y, d_moj.(_ar(rg))') # outer-product
end

_tau_arc(dsd::RealU, ds::RealU, tangam::Real) =
    dsd / ds * atan(tangam)
_tau_flat(dsd::RealU, ds::RealU, tangam::Real) =
    dsd / ds * tangam
_tau(rg::Union{SinoFanArc,CtFanArc}, tangam) =
    _tau_arc.(rg.dsd, rg.d, tangam)
_tau(rg::Union{SinoFanFlat,CtFanFlat}, tangam) =
    _tau_flat.(rg.dsd, rg.d, tangam)

function _tau(
    rg::Union{SinoFan,CtFan},
    x::AbstractVector,
    y::AbstractVector,
)
    b = _ar(rg)' # row vector, for outer-product
    xb = x * cos.(b) + y * sin.(b)
    yb = -x * sin.(b) + y * cos.(b)
    if rg isa SinoFan
        xb .-= rg.source_offset
    end
    tangam = xb ./ (_dso(rg) .- yb) # e,tomo,fan,L,gam
    tau = _tau(rg, tangam)
    return tau
end


"""
    _tau(rg::RayGeom, x, y)
Projected `s/ds`, useful for footprint center and support.
Returns `Matrix` of size `length(x) × rg.na`.
"""
function _tau(
    rg::RayGeom,
    x::AbstractArray,
    y::AbstractArray,
)
    size(x) == size(y) || throw("bad x,y size")
    return _tau(rg, vec(x), vec(y))
end

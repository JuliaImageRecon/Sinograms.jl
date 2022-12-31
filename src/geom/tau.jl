#=
geom/tau.jl
Projection position of pixel centers in parallel-beam coordinates.
=#


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
    return _tau.(_ar(rg)', x, y, _ds(rg)) # outer-product
end

# This one may not be useful but it helps unify tests:
function _tau(rg::SinoMoj, x::AbstractVector, y::AbstractVector)
    d_moj = _d_moj(rg)
    return _tau.(_ar(rg)', x, y, d_moj.(_ar(rg))') # outer-product
end

_tau_arc(dsd::RealU, ds::RealU, tangam::Real) =
    dsd / ds * atan(tangam)

_tau_flat(dsd::RealU, ds::RealU, tangam::Real) =
    dsd / ds * tangam

_tau(rg::Union{SinoFanArc,CtFanArc}, tangam) =
    _tau_arc.(rg.dsd, _ds(rg), tangam)

_tau(rg::Union{SinoFanFlat,CtFanFlat}, tangam) =
    _tau_flat.(rg.dsd, _ds(rg), tangam)

function _tau(
    rg::Union{SinoFan,CtFan},
    x::AbstractVector,
    y::AbstractVector,
)
    axes(x) == axes(y) || error("x,y axes")
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
    axes(x) == axes(y) || error("x,y axes")
    return _tau(rg, vec(x), vec(y))
end


_tau(rg::RayGeom) = (x,y) -> _tau(rg, x, y)

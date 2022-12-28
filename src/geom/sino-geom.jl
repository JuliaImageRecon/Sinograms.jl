#=
geom/sino-geom.jl
sinogram geometry definitions for 2D tomographic image reconstruction
2019-07-01, Jeff Fessler, University of Michigan
2022-01-22, copied from MIRT.jl and updated to support Unitful values
=#

export dims, downsample, oversample, rays, axes
export sino_w, sino_s

# Methods common to all types

dims(sg::SinoGeom) = (sg.nb, sg.na)::NTuple{2,Int}

# use ° via `angles()` because mainly for plots
Base.axes(sg::SinoGeom) = (sino_s(sg), angles(sg))

sino_w(sg::SinoGeom) = Toffset((sg.nb-1)/2 + sg.offset)::Toffset

#_geom_s(sg::SinoGeom) = sg.d * ((0:sg.nb-1) .- sino_w(sg))
#sino_s(sg::SinoGeom) = sg.d * ((0:sg.nb-1) .- sino_w(sg))
#sino_s(sg::SinoGeom) = _lin_range(sg.d, sg.w, sg.nb) # can't infer!?
function sino_s(
    st::SinoGeom{Td} ;
    T::Type{<:Number} = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.d, st.w, st.nb)
end
_geom_s(st::SinoGeom) = sino_s(st) # todo
_s(st::SinoGeom) = sino_s(st) # todo

sino_dr(sg::SinoGeom) = sg.d
sino_dr(sg::SinoMoj) = NaN



# down/up sampling

function _downsample(sg::SinoGeom, down::Int)
    nb = 2 * max(sg.nb ÷ 2down, 1) # keep it even
    na = max(sg.na ÷ down, 1)
    return (nb, sg.d * down, sg.offset, na, sg.orbit, sg.orbit_start)
end


"""
    downsample(st, down)

Down-sample sinogram geometry
(for testing with small problems).
"""
function downsample(st::G, down::Int) where {G <: SinoParallel}
    return (down == 1) ? st : G(_downsample(st, down)...)::G
end

function downsample(st::G, down::Int) where {G <: SinoFan}
    return (down == 1) ? st :
        G(_downsample(st, down)..., st.source_offset, st.dsd, st.dod)::G
end


function _oversample(st::SinoGeom, over::Int)
    return (
        st.nb * over, st.d / over, st.offset * over,
        st.na, st.orbit, st.orbit_start,
    )
end

"""
    oversample(st, over::Int)

Over-sample sinogram geometry in "radial" dimension.
For Mojette sampling, it means that `d = dx/over`.
"""
function oversample(st::G, over::Int) where {G <: SinoParallel}
    return (over == 1) ? st : G(_oversample(st, over)...)::G
end

function oversample(st::G, over::Int) where {G <: SinoFan}
    return (over == 1) ? st :
        G(_oversample(st, over)..., st.source_offset, st.dsd, st.dod)::G
end


#=
# gamma for general finite dfs (unsupported, would need new SinoFan subtype with dfs)
function sino_geom_gamma(sg::SinoFan)
    dis_foc_det = sg.dfs + sg.dsd
    α = sino_s(sg) / dis_foc_det # equivalent to s/dsd when dfs=0
    return atan.(dis_foc_det * sin.(α), dis_foc_det * cos.(α) .- sg.dfs)
end
=#




# type inference help:
function _rays_type2(Td,To)
    Tϕ = eltype(oneunit(to_radians([oneunit(To)])[1]))
    return Iterators.ProductIterator{Tuple{
        LinRange{eltype(1f0 * oneunit(Td)), Int},
        LinRange{Tϕ, Int},
    }}
end


"""
     = rays(sg::SinoGeom)

Radial `r` and angular `ϕ` coordinates (in radians)
of all sinogram elements
for the given geometry.
Return type of `i` is a `ProductIterator` that makes tuples of the form
`(r, ϕ)`.
To make projections call
`p = [fun(c...) for c in i]` where `fun` is `radon(...)`.
"""
function rays(st::SinoPar{Td,To})::_rays_type2(Td,To) where {Td,To}
    s = sino_s(st)
    ϕ = st.ar # / oneunit(eltype(st.ar)) # deg2rad.(angles(st))
    i = Iterators.product(s, ϕ)
    return i
end


_moj_to_par(s, ϕ) = (s * maximum(abs, sincos(ϕ)), ϕ)
_moj_to_par(sϕ::Tuple) = _moj_to_par(sϕ...)
_moj_rays_type2(Td,To) = Base.Generator{_rays_type2(Td,To), typeof(_moj_to_par)}

function rays(st::SinoMoj{Td,To})::_moj_rays_type2(Td,To) where {Td, To}
    s = sino_s(st)
    ϕ = st.ar # / oneunit(eltype(st.ar))
    i = Iterators.product(s, ϕ)
    return Iterators.map(_moj_to_par, i)
end


function _fan_to_par(st::SinoFan, s, β)
    γ = _geom_gamma_s(st, s)
    r = st.dso * sin(γ)
    ϕ = γ + β
    return (r, ϕ)
end
_fan_to_par(st::SinoFan, sβ::Tuple) = _fan_to_par(st, sβ...)

function rays(st::SinoFan{Td,To}) where {Td, To}
    s = sino_s(st)
    β = st.ar # / oneunit(eltype(st.ar))
    i = Iterators.product(s, β)
    fun = sβ -> _fan_to_par(st, sβ) # closure prevents type inference?
    return Iterators.map(fun, i)
end



# τ (unitless)
function _sino_geom_tau(ϕ::RealU, x::RealU, y::RealU, dr::RealU)
    sϕ, cϕ = sincos(ϕ)
    return (x * cϕ + y * sϕ) / dr
end

function _sino_geom_tau(
    st::Union{SinoPar,CtPar},
    x::AbstractVector,
    y::AbstractVector,
)
    return _sino_geom_tau.(st.ar', x, y, st.dr) # outer-product
end

# this one may not be useful but it helps unify tests
_sino_geom_tau(sg::SinoMoj, x::AbstractVector, y::AbstractVector) =
    _sino_geom_tau.(sg.ar', x, y, sg.d_moj.(sg.ar)') # outer-product

_sino_geom_tau_arc(dsd::RealU, ds::RealU, tangam::Real) =
    dsd / ds * atan(tangam)
_sino_geom_tau_flat(dsd::RealU, ds::RealU, tangam::Real) =
    dsd / ds * tangam
_sino_geom_tau(st::Union{SinoFanArc,CtFanArc}, tangam) =
    _sino_geom_tau_arc.(st.dsd, st.ds, tangam)
_sino_geom_tau(st::Union{SinoFanFlat,CtFanFlat}, tangam) =
    _sino_geom_tau_flat.(st.dsd, st.ds, tangam)

function _sino_geom_tau(
    st::Union{SinoFan,CtFan},
    x::AbstractVector,
    y::AbstractVector,
)
    b = st.ar' # row vector, for outer-product
    xb = x * cos.(b) + y * sin.(b)
    yb = -x * sin.(b) + y * cos.(b)
    if st isa SinoFan
        xb .-= st.source_offset
    end
    tangam = xb ./ (st.dso .- yb) # e,tomo,fan,L,gam
    tau = _sino_geom_tau(st, tangam)
    return tau
end


"""
    sino_geom_tau(st::Union{SinoGeom,CtGeom}, x, y)
Projected `s/ds`, useful for footprint center and support.
Returns `Matrix` of size `length(x) × st.na`.
"""
function sino_geom_tau(
    st::Union{SinoGeom{Td,To},CtGeom{Td,To}},
    x::AbstractArray{Tx},
    y::AbstractArray{Ty},
) where {Td <: Number, To <: Number, Tx <: Number, Ty <: Number}
    size(x) == size(y) || throw("bad x,y size")
    T = promote_type(map(t -> eltype(one(t)), (Tx, Ty, Td, To))...)
    return _sino_geom_tau(st, vec(x), vec(y))::Matrix{T}
end

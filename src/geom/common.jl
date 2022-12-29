#=
geom/common.jl
Methods that are common to both 2d and 3d
=#

export angles


"""
    show(io::IO, ::MIME"text/plain", rg::RayGeom)
"""
Base.show(io::IO, ::MIME"text/plain", rg::RayGeom) = _show(io, MIME("text/plain"), rg)

Base.ones(T::Type{<:Number}, rg::RayGeom) = ones(T, dims(rg))
Base.ones(rg::RayGeom) = ones(Float32, rg)
Base.zeros(T::Type{<:Number}, rg::RayGeom) = zeros(T, dims(rg))
Base.zeros(rg::RayGeom) = zeros(Float32, rg)

"""
    angles(rg::RayGeom) =
Return vector of angles
for this ray geometry,
in whatever units the user used to specify
`orbit` and `orbit_start`,
typically degrees.
"""
angles(rg::RayGeom) =
    range(rg.orbit_start, length = rg.na, step = rg.orbit / rg.na)

# angles in radians
_ar(rg::RayGeom) = to_radians(angles(rg))

# minimum orbit for a fan-beam short scan
_orbit_short(rg::RayGeom) = 180 + 2 * rad2deg(_gamma_max(rg)) # (degrees)

# distance from detector arc focal spot to source for fan-beam
_dfs(::Union{SinoFanArc{Td},CtFanArc{Td}}) where Td = zero(Td)
_dfs(::Union{SinoFanFlat{Td},CtFanFlat{Td}}) where Td = Inf * oneunit(Td)

# distance from source to origin for fan-beam
_dso(rg::Union{SinoFan,CtFan}) = rg.dsd - rg.dod

"""
    _rfov(rg::RayGeom)
Radial FOV.
"""
_rfov(rg::Union{SinoPar,CtPar}) = maximum(abs, _s(rg))
_rfov(rg::SinoMoj) = dims(rg)[1]/2 * minimum(_d_ang(rg)) # (ignores offset)
_rfov(rg::Union{SinoFan,CtFan}) = _dso(rg) * sin(_gamma_max(rg))


"""
    _xds(rg::RayGeom)
Center `x` positions of detectors (for beta = 0),
for central row of detector.
"""
_xds(rg::Union{SinoPar,CtPar}) = _s(rg)
_xds(rg::SinoMoj) = _s(rg) # todo: really should be angle dependent
_xds(rg::Union{SinoFanArc,CtFanArc}) = rg.dsd * sin.(_gamma(rg)) .+ rg.source_offset
_xds(rg::Union{SinoFanFlat,CtFanFlat}) = _s(rg) .+ rg.source_offset


"""
    _yds(rg::RayGeom)
Center `y` positions of detectors (for beta = 0),
for central row of detector.
"""
_yds(rg::Union{SinoPar{Td},CtPar{Td}}) where Td = zeros(Td, dims(rg)[1])
_yds(rg::SinoMoj{Td}) where Td = zeros(Td, dims(rg)[1])
_yds(rg::Union{SinoFanArc,CtFanArc}) = _dso(rg) .- rg.dsd * cos.(_gamma(rg))
_yds(rg::Union{SinoFanFlat,CtFanFlat}) = fill(-rg.dod, dims(rg)[1])



"""
    _gamma(rg::RayGeom [, s])
Return gamma (γ: fan angle) values for fan-beam geometry.
"""
_gamma

_gamma(rg::Union{SinoFanArc,CtFanArc}, s::RealU) = s / rg.dsd
_gamma(rg::Union{SinoFanFlat,CtFanFlat}, s::RealU) = atan(s / rg.dsd)
_gamma(rg::Union{SinoFanArc,CtFanArc}, ss::AbstractArray) = ss / rg.dsd
_gamma(rg::Union{SinoFanFlat,CtFanFlat}, ss::AbstractArray) = @. atan(ss / rg.dsd)
_gamma(rg::Union{SinoFan, CtFan}) = _gamma(rg, _s(rg))

_gamma_max(rg::RayGeom) = maximum(_gamma(rg))
_gamma_max_abs(rg::RayGeom) = maximum(abs, _gamma(rg))


"""
    _shape(rg, x::AbstractArray)
Reshape `x` to `dims(rg)` or `(dims(rg)..., :)`.
Not type stable.
"""
_shape(rg::RayGeom, x::AbstractArray) = _shape(x, dims(rg))


"""
    _unitv([T=Float32], rg:RayGeom [, pos::Tuple])
Projection views with a single non-zero ray value
at position `pos` (default: middle).
"""
function _unitv(
    T::Type{<:Number},
    rg::RayGeom,
    pos::Tuple = dims(rg) .÷ 2 .+ 1,
)
    out = zeros(T, rg)
    out[pos...] = one(T)
    return out
end

_unitv(rg::RayGeom, args... ; kwargs...) = _unitv(Float32, rg, args... ; kwargs...)


_taufun(rg::RayGeom) = (x,y) -> _tau(rg, x, y)

_d_moj(rg::SinoMoj) = ar -> rg.d * max(abs(cos(ar)), abs(sin(ar)))
_d_ang(rg::SinoMoj) = _d_moj(rg).(_ar(rg))

#=
geom/common.jl
Methods that are common to both 2d and 3d
=#

export angles


"""
    show(io::IO, ::MIME"text/plain", st::RayGeom)
"""
Base.show(io::IO, ::MIME"text/plain", st::RayGeom) = _show(io, MIME("text/plain"), st)


# todo: for now just SinoGeom
function _getproperty(st::SinoGeom, s::Symbol, arg) # not type stable :(
    d = Dict(arg)
    return haskey(d, s) ? d[s](st) : getfield(st, s)
end

Base.getproperty(st::SinoGeom, s::Symbol) = _getproperty(st, s, _props(st))

Base.propertynames(st::SinoGeom) =
    (fieldnames(typeof(st))..., map(x -> x[1], _props(st))...)

#=
function _getproperty(st::RayGeom, s::Symbol, arg) # not type stable :(
    d = Dict(arg)
    return haskey(d, s) ? d[s](st) : getfield(st, s)
end

Base.getproperty(st::RayGeom, s::Symbol) = _getproperty(st, s, _props(st))

Base.propertynames(st::RayGeom) =
    (fieldnames(typeof(st))..., map(x -> x[1], _props(st))...)
=#


Base.ones(T::Type{<:Number}, st::RayGeom) = ones(T, dims(st))
Base.ones(st::RayGeom) = ones(Float32, st)
Base.zeros(T::Type{<:Number}, st::RayGeom) = zeros(T, dims(st))
Base.zeros(st::RayGeom) = zeros(Float32, st)

#angles(st::RayGeom{Td,To}) where {Td,To} =
angles(st::RayGeom) =
    range(st.orbit_start, length = st.na, step = st.orbit / st.na)
_ar(st::RayGeom) = to_radians(angles(st))
# todo type inference issues
#   LinRange(st.orbit_start, To(st.orbit_start + (st.na-1)/st.na * st.orbit), st.na)::LinRange{To,Int}

_orbit_short(st::RayGeom) = 180 + 2 * rad2deg(_gamma_max(st)) # (degrees)

_dfs(st::Union{SinoFanArc{Td},CtFanArc{Td}}) where Td = zero(Td)
_dfs(st::Union{SinoFanFlat{Td},CtFanFlat{Td}}) where Td = Inf * oneunit(Td)
_geom_dfs = _dfs # todo cut

_dso(st::Union{SinoFan,CtFan}) = st.dsd - st.dod

"""
    _rfov(st::RayGeom)
Radial FOV.
"""
_rfov(st::Union{SinoPar,CtPar}) = maximum(abs, _s(st))
_rfov(st::SinoMoj) = dims(st)[1]/2 * minimum(st.d_ang) # (ignores offset)
_rfov(st::Union{SinoFan,CtFan}) = _dso(st) * sin(_gamma_max(st))
_geom_rfov = _rfov # todo cut


"""
    _xds(st::RayGeom)
Center `x` positions of detectors (for beta = 0),
for central row of detector.
"""
_xds(st::Union{SinoPar,CtPar}) = _s(st)
_xds(st::SinoMoj) = _s(st) # todo: really should be angle dependent
_xds(st::Union{SinoFanArc,CtFanArc}) = st.dsd * sin.(_gamma(st)) .+ st.source_offset
_xds(st::Union{SinoFanFlat,CtFanFlat}) = _s(st) .+ st.source_offset
_geom_xds = _xds # todo cut


"""
    _yds(st::RayGeom)
Center `y` positions of detectors (for beta = 0),
for central row of detector.
"""
_yds(st::Union{SinoPar{Td},CtPar{Td}}) where Td = zeros(Td, dims(st)[1])
_yds(st::SinoMoj{Td}) where Td = zeros(Td, dims(st)[1])
_yds(st::Union{SinoFanArc,CtFanArc}) = _dso(st) .- st.dsd * cos.(_gamma(st))
_yds(st::Union{SinoFanFlat,CtFanFlat}) = fill(-st.dod, dims(st)[1])
_geom_yds = _yds # todo cut



"""
    _gamma(st::RayGeom [, s])
Return gamma (γ: fan angle) values for fan-beam geometry.
"""
_gamma

_gamma(st::Union{SinoFanArc,CtFanArc}, s::RealU) = s / st.dsd
_gamma(st::Union{SinoFanFlat,CtFanFlat}, s::RealU) = atan(s / st.dsd)
_gamma(st::Union{SinoFanArc,CtFanArc}, ss::AbstractArray) = ss / st.dsd
_gamma(st::Union{SinoFanFlat,CtFanFlat}, ss::AbstractArray) = @. atan(ss / st.dsd)
_gamma(st::Union{SinoFan, CtFan}) = _gamma(st, _s(st))
_geom_gamma = _gamma # todo cut

#=
_gamma(st::Union{SinoFanArc, CtFanArc}) = _s(st) / st.dsd
_gamma(st::Union{SinoFanFlat, CtFanFlat}) = atan.(_s(st) / st.dsd)

function _gamma(st::Union{SinoFanArc{Td}, CtFanArc{Td}}) where Td # todo cut Td
#   s = sino_s(st)
    s = _geom_s(st)
#   s = st.s
    gamma = s / st.dsd # 3rd gen: equiangular
#   Tg = eltype(one(Td))
    return gamma#::LinRange{Tg,Int}
end

function _gamma(st::Union{SinoFanFlat{Td}, CtFanFlat{Td}}) where Td
#   s = sino_s(st)
#   s = st.s
    s = _geom_s(st)
    gamma = atan.(s / st.dsd) # flat
#   Tg = eltype(one(Td))
    return gamma#::Vector{Tg}
end

_geom_gamma_s(st::Union{SinoFanArc,CtFanArc}, s::RealU) = s / st.dsd
_geom_gamma_s(st::Union{SinoFanFlat,CtFanFlat}, s::RealU) = atan(s / st.dsd)
_geom_gamma_s(st::Union{SinoFanArc,CtFanArc}, ss::AbstractArray) = ss / st.dsd
_geom_gamma_s(st::Union{SinoFanFlat,CtFanFlat}, ss::AbstractArray) = atan.(ss / st.dsd)
#_geom_gamma_s(st::RayGeom, ss::AbstractArray) = (s::RealU -> _geom_gamma_s(st, s)).(ss)
#_geom_gamma(st::RayGeom) = _geom_gamma_s(st, st.s)
=#
_geom_gamma_s = _gamma # todo cut

_gamma_max(st::RayGeom) = maximum(_gamma(st))
_gamma_max_abs(st::RayGeom) = maximum(abs, _gamma(st))
#_geom_gamma_max = _gamma_max # todo cut
#_geom_gamma_max_abs = _geom_gamma_max_abs # todo cut

"""
    _shape(st, x::AbstractArray)
Reshape `x` to `dims(st)` or `(dims(st)..., :)`.
"""
_shape(st::RayGeom, x::AbstractArray) = _shape(x, dims(st))


"""
    _unitv([T=Float32], st:RayGeom [, pos::Tuple])
Projection views with a single non-zero ray value
at position `pos` (default: middle).
"""
function _unitv(
    T::Type{<:Number},
    st::RayGeom,
    pos::Tuple = dims(st) .÷ 2 .+ 1,
)
    out = zeros(T, st)
    out[pos...] = one(T)
    return out
end

_unitv(st::RayGeom, args... ; kwargs...) = _unitv(Float32, st, args... ; kwargs...)
_geom_unitv = _unitv # todo cut

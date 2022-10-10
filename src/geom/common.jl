#=
geom/common.jl
Methods that are common to both 2d and 3d
=#

export angles


"""
    show(io::IO, ::MIME"text/plain", st::RayGeom)
"""
Base.show(io::IO, ::MIME"text/plain", st::RayGeom) = _show(io,  MIME("text/plain"), st)


function _getproperty(st::RayGeom, s::Symbol, arg)
    d = Dict(arg)
    return haskey(d, s) ? d[s](st) : getfield(st, s)
end

Base.getproperty(st::RayGeom, s::Symbol) = _getproperty(st, s, _props(st))

Base.propertynames(st::RayGeom) =
    (fieldnames(typeof(st))..., map(x -> x[1], _props(st))...)


Base.ones(T::DataType, st::RayGeom) = ones(T, dims(st))
Base.ones(st::RayGeom) = ones(Float32, st)
Base.zeros(T::DataType, st::RayGeom) = zeros(T, dims(st))
Base.zeros(st::RayGeom) = zeros(Float32, st)

angles(st::RayGeom{Td,To}) where {Td,To} =
    LinRange(st.orbit_start, To(st.orbit_start + (st.na-1)/st.na * st.orbit), st.na)::LinRange{To,Int}
#   range(st.orbit_start, length = st.na, step = st.orbit / st.na) # type inference issues

_orbit_short(st::RayGeom) = 180 + 2 * rad2deg(st.gamma_max) # (degrees)

_geom_dfs(st::Union{SinoFanArc{Td},CtFanArc{Td}}) where Td = zero(Td)
_geom_dfs(st::Union{SinoFanFlat{Td},CtFanFlat{Td}}) where Td = Inf * oneunit(Td)


"""
    _geom_rfov(st::RayGeom)
Radial FOV.
"""
_geom_rfov(st::Union{SinoPar,CtPar}) = maximum(abs, st.s)
_geom_rfov(st::SinoMoj) = dims(st)[1]/2 * minimum(st.d_ang) # (ignores offset)
_geom_rfov(st::Union{SinoFan,CtFan}) = st.dso * sin(st.gamma_max)


"""
    _geom_xds(st::RayGeom)
Center `x` positions of detectors (for beta = 0),
for central row of detector.
"""
_geom_xds(st::Union{SinoPar,CtPar}) = st.s
_geom_xds(st::SinoMoj) = st.s # todo: really should be angle dependent
_geom_xds(st::Union{SinoFanArc,CtFanArc}) = st.dsd * sin.(st.gamma) .+ st.source_offset
_geom_xds(st::Union{SinoFanFlat,CtFanFlat}) = st.s .+ st.source_offset


"""
    _geom_yds(st::RayGeom)
Center `y` positions of detectors (for beta = 0),
for central row of detector.
"""
_geom_yds(st::Union{SinoPar{Td},CtPar{Td}}) where Td = zeros(Td, dims(st)[1])
_geom_yds(st::SinoMoj{Td}) where Td = zeros(Td, dims(st)[1])
_geom_yds(st::Union{SinoFanArc,CtFanArc}) = st.dso .- st.dsd * cos.(st.gamma)
_geom_yds(st::Union{SinoFanFlat,CtFanFlat}) = fill(-st.dod, dims(st)[1])



"""
    _geom_gamma(st::RayGeom)
Return gamma (γ: fan angle) values for fan-beam geometry.
"""
function _geom_gamma(st::Union{SinoFanArc{Td}, CtFanArc{Td}}) where Td
#   s = sino_s(st)
    s = st.s
    gamma = s / st.dsd # 3rd gen: equiangular
    Tg = eltype(one(Td))
    return gamma::LinRange{Tg,Int}
end

function _geom_gamma(st::Union{SinoFanFlat{Td}, CtFanFlat{Td}}) where Td
#   s = sino_s(st)
    s = st.s
    gamma = atan.(s / st.dsd) # flat
    Tg = eltype(one(Td))
    return gamma::Vector{Tg}
end

_geom_gamma_s(st::Union{SinoFanArc,CtFanArc}, s::RealU) = s / st.dsd
_geom_gamma_s(st::Union{SinoFanFlat,CtFanFlat}, s::RealU) = atan(s / st.dsd)
_geom_gamma_s(st::Union{SinoFanArc,CtFanArc}, ss::AbstractArray) = ss / st.dsd
_geom_gamma_s(st::Union{SinoFanFlat,CtFanFlat}, ss::AbstractArray) = atan.(ss / st.dsd)
#_geom_gamma_s(st::RayGeom, ss::AbstractArray) = (s::RealU -> _geom_gamma_s(st, s)).(ss)

#_geom_gamma(st::RayGeom) = _geom_gamma_s(st, st.s)
_geom_gamma_max(st::RayGeom) = maximum(_geom_gamma(st))
_geom_gamma_max_abs(st::RayGeom) = maximum(abs, _geom_gamma(st))


"""
    _geom_unitv([T=Float32], st:RayGeom [, pos::Tuple])
Projection views with a single non-zero ray value
at position `pos` (default: middle).
"""
function _geom_unitv(
    T::DataType,
    st::RayGeom,
    pos::Tuple = dims(st) .÷ 2 .+ 1,
)
    out = zeros(T, st)
    out[pos...] = one(T)
    return out
end

_geom_unitv(st::RayGeom, args... ; kwargs...) = _geom_unitv(Float32, st, args... ; kwargs...)

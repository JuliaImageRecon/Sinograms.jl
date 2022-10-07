# fbp/ramp.jl

export fbp_ramp, ramp_flat, ramp_arc

#using Sinograms: SinoFanFlat, SinoFanArc, SinoPar, CtGeom, RealU


"""
    h, n = fbp_ramp(sg::SinoGeom, N::Int)
'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
This sampled band-limited approach avoids the aliasing that would be
caused by sampling the ramp directly in the frequency domain.

in
- `sg::SinoGeom`
- `N::Int` : # of samples (must be even)

out
- `h::Vector{<:RealU}` : samples of band-limited ramp filter
- `n::UnitRange{Int64}` : -(N÷2):(N÷2-1)
"""
function fbp_ramp(sg::SinoPar{Td}, N::Int) where Td
    R = _ramp_type(Td)
    T = Tuple{Vector{R}, UnitRange{Int64}}
    return ramp_flat(N, sg.d)::T
end

# for Mojette, `dr` varies with projection angle
function fbp_ramp(sg::SinoMoj{Td}, N::Int ; dr::Td = sg.d) where Td
    R = _ramp_type(Td)
    T = Tuple{Vector{R}, UnitRange{Int64}}
    return ramp_flat(N, dr)::T
end

function fbp_ramp(sg::Union{SinoFanFlat{Td},CtFanFlat{Td}}, N::Int) where Td
    R = _ramp_type(Td)
    T = Tuple{Vector{R}, UnitRange{Int64}}
    return ramp_flat(N, sg isa SinoGeom ? sg.d : sg.ds)::T
end

function fbp_ramp(sg::Union{SinoFanArc{Td},CtFanArc{Td}}, N::Int) where Td
    R = _ramp_type(Td)
    T = Tuple{Vector{R}, UnitRange{Int64}}
   return ramp_arc(N, sg isa SinoGeom ? sg.d : sg.ds, sg.dsd)::T
end

function _ramp_type(arg...)
    R = promote_type(arg...)
    R = promote_type(R, eltype(1f0 * oneunit(R))) # at least Float32
    R = eltype(1 / oneunit(R)^2)
    return R
end

function _ramp_arc(n::Int, ds::RealU, dsd::RealU)
    R = _ramp_type(eltype(ds), eltype(dsd))
    h = n == 0 ? 0.25 / abs2(ds) :
        isodd(n) ? -1 / abs2(π * dsd * sin(n * ds / dsd)) :
        zero(1 / abs2(ds))
    return R(h)
end


function _ramp_flat(n::Int, ds::RealU)
    R = _ramp_type(eltype(ds))
    h = n == 0 ? 0.25 / abs2(ds) :
        isodd(n) ? -1 / abs2(π * n * ds) :
        zero(1 / abs2(ds))
    return R(h)
end


function ramp_arc(N::Int, ds::RealU, dsd::RealU)
    isodd(N) && throw("N must be even")

    (N/2 * ds / dsd > 0.9 * π/2) &&
        throw("angle is $(rad2deg(N/2 * ds / dsd)) degrees: too large;
        physically impossible arc geometry")

    n = -(N÷2):(N÷2-1)
    h = _ramp_arc.(n, ds, dsd)
    return h, n
end


function ramp_flat(N::Int, ds::RealU)
    isodd(N) && throw("N must be even")
    n = -(N÷2):(N÷2-1)
    h = _ramp_flat.(n, ds)
    return h, n
end

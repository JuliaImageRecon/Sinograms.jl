# fbp/ramp.jl

export fbp_ramp, ramp_flat, ramp_arc

#using Sinograms: SinoFanFlat, SinoFanArc, SinoPar, CtGeom, RealU


"""
    h, n = fbp_ramp(rg::SinoGeom, N::Int)
'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
This sampled band-limited approach avoids the aliasing that would be
caused by sampling the ramp directly in the frequency domain.

in
- `rg::SinoGeom`
- `N::Int` : # of samples (must be even)

out
- `h::Vector{<:RealU}` : samples of band-limited ramp filter
- `n::UnitRange{Int64}` : -(N÷2):(N÷2-1)
"""
fbp_ramp

fbp_ramp(rg::SinoPar, N::Int) = ramp_flat(N, rg.d)

# for Mojette, `dr` varies with projection angle
function fbp_ramp(rg::SinoMoj{Td}, N::Int ; dr::Td = rg.d) where Td
    return ramp_flat(N, dr)
end

fbp_ramp(rg::Union{SinoFanFlat,CtFanFlat}, N::Int) =
    ramp_flat(N, _ds(rg))

fbp_ramp(rg::Union{SinoFanArc,CtFanArc}, N::Int) =
    ramp_arc(N, _ds(rg), rg.dsd)

function _ramp_type(arg...)
    R = promote_type(arg...)
    R = promote_type(R, typeof(1f0 * oneunit(R))) # at least Float32
    R = eltype(1 / oneunit(R)^2)
    return R
end

function _ramp_arc(n::Int, ds::RealU, dsd::RealU)
    R = _ramp_type(eltype(ds), eltype(dsd))
    h = n == 0 ? R(0.25 / abs2(ds)) :
        isodd(n) ? R(-1 / abs2(π * dsd * sin(n * ds / dsd))) :
        zero(R)
    return R(h)
end

function _ramp_flat(n::Int, ds::RealU)
    R = _ramp_type(eltype(ds))
    h = n == 0 ? R(0.25 / abs2(ds)) :
        isodd(n) ? R(-1 / abs2(π * n * ds)) :
        zero(R)
    return h
end


"""
    (h, n) = ramp_arc(N::Int, ds::RealU, dsd::RealU)
Ramp filter samples for arc fan geometry,
for `n = -(N÷2):(N÷2-1)`.
- `N` must be even.
"""
function ramp_arc(N::Int, ds::RealU, dsd::RealU)
    isodd(N) && throw("N must be even")

    (N/2 * ds / dsd > 0.9 * π/2) &&
        throw("angle is $(rad2deg(N/2 * ds / dsd)) degrees: too large;
        physically impossible arc geometry")

    n = -(N÷2):(N÷2-1)
    h = _ramp_arc.(n, ds, dsd)
    return h, n
end


"""
    (h, n) = ramp_flat(N::Int, ds::RealU)
Ramp filter samples for flat fan geometry,
for `n = -(N÷2):(N÷2-1)`.
- `N` must be even.
"""
function ramp_flat(N::Int, ds::RealU)
    isodd(N) && throw("N must be even")
    n = -(N÷2):(N÷2-1)
    h = _ramp_flat.(n, ds)
    return h, n
end

#=
fbp/parker.jl
Parker weighting for FBP reconstruction
http://doi.org/10.1118/1.595078
=#

export parker_weight
# using Sinograms: SinoGeom, SinoPar, SinoFan, SinoMoj
using LazyGrids: ndgrid


"""
    parker_weight(rg::SinoGeom, T = Float32)
Compute Parker weighting for non-360° orbits.
See http://doi.org/10.1118/1.595078.
Returns `Matrix{T}` of size:
- (1,1) for `SinoPar` with typical 180 or 360 orbit
- (1,na) for `SinoPar` with atypical orbit
- (1,1) for `SinoFan` with typical 360 orbit
- (ns,na) for `SinoFan` with typical 360 orbit
"""
parker_weight


# parallel-beam case
function parker_weight_par(
    orbit::RealU,
    ad::AbstractVector{<:RealU}, # angles in degrees: angles(rg) - rg.orbit_start
    T::Type{<:AbstractFloat} = Float32,
)

    da = diff(ad)
    all(≈(da[begin]), da) || error("Unequal angles not yet supported")

    if (orbit ÷ 180) * 180 == orbit
        return ones(T, 1, 1) # no weighting needed
    end

    if abs(orbit) < 180
        @warn("orbit $orbit < 180")
        return ones(T, 1, 1) # nonuniform weighting would not help
    end

    orbit = abs(orbit)
    orbit > 360 && error("only 180 ≤ |orbit| ≤ 360 supported for Parker weighting")
    extra = orbit - 180 # extra beyond 180

    wt = ones(T, 1, length(ad)) # (1,na)
    ad = abs.(ad)
    ii = ad .< extra
    @. wt[ii] = abs2(sin(ad[ii] / extra * π/2))
    ii = ad .≥ 180
    @. wt[ii] = abs2(sin((orbit - ad[ii]) / extra * π/2))
    return wt
end


parker_weight(rg::SinoPar, T::Type{<:AbstractFloat} = Float32) =
    parker_weight_par(rg.orbit, angles(rg) .- rg.orbit_start, T)


function parker_weight_fan_short(
    nb::Int,
    na::Int,
    orbit::RealU,
    orbit_short::RealU,
    ar::AbstractVector{<:RealU}, # angles in radians
    gam::AbstractVector{<:RealU},
    gamma_max::RealU, # half of fan angle
    T::Type{<:AbstractFloat} = Float32,
)

    da = diff(ar)
    all(≈(da[begin]), da) || error("Unequal angles not yet supported")
#   Base.require_one_based_indexing(ar, gam)

    na == length(ar) || error("na bug")
    orbit < orbit_short - eps() &&
        @warn("orbit $orbit is less than a short scan $orbit_short")

    orbit > orbit_short + rad2deg(ar[begin+1] - ar[begin]) &&
        @warn("orbit $orbit exceeds short scan $orbit_short by %g views")
    #   (orbit - orbit_short) / rad2deg(ar[begin+1] - ar[begin]))

    bet = ar .- ar[begin] # trick: force 0 start, so this ignores orbit_start!
    bet = max.(bet, zero(eltype(bet))) # fix for issue #49
    (gg, bb) = ndgrid(gam, bet)

    fun = (x) -> sin(π/2 * x)^2 # smooth out [0,1] ramp
    #=
    todo: We could use integrals of this function
    over the tiny angular range of each projection view,
    so that the sum over beta of these functions is a constant.
    =#

    wt = zeros(T, nb, na) # any extra views will get 0 weight
    ii = @. bb < 2 * (gamma_max .- gg) # 0 <= bb not needed
    tmp = bb[ii] ./ (2 * (gamma_max .- gg[ii]))
    @. wt[ii] = fun(tmp)

    ii = @. (2 * (gamma_max - gg) ≤ bb) & (bb < π - 2 * gg)
    wt[ii] .= 1

    ii = @. (π - 2 * gg < bb) & (bb ≤ π + 2 * gamma_max)
    tmp = @. (π + 2*gamma_max - bb[ii]) / (2 * (gamma_max + gg[ii]))
    @. wt[ii] = fun(tmp)
    return wt
end


function parker_weight(rg::SinoFan, T::Type{<:AbstractFloat} = Float32)
    if (rg.orbit ÷ 360) * 360 == rg.orbit
        return ones(T, 1, 1) # no weighting needed
    end
    return parker_weight_fan_short(
        rg.nb, rg.na, rg.orbit, _orbit_short(rg),
        _ar(rg), _gamma(rg), _gamma_max(rg) #; T
    )
end


function parker_weight(rg::SinoMoj ; T::Type{<:AbstractFloat} = Float32)
    orbit = abs(rg.orbit)
    na = rg.na
    ((rg.orbit ÷ 180) * 180 == orbit) ||
        throw("No Parker weighting for Mojette geometry with orbit=$orbit")
    wt = ones(T, 1, 1)
    return wt
end


function parker_weight_fan_short(rg::CtFan; kwargs...)
    weight = parker_weight_fan_short(
        rg.ns, rg.na, rg.orbit, _orbit_short(rg),
        _ar(rg), _gamma(rg), _gamma_max(rg); kwargs...,
    )
    return weight
end


"""
    parker_weight(rg::CtFan; T::Type{<:AbstractFloat} = Float32, kwargs...)
For 3D case, return `Array{T,3}` where size is
- `(1,1,1)` typical fan case with 360° orbit
- `(ns,1,na)` atypical fan case including short scan
"""
function parker_weight(rg::CtFan; T::Type{<:AbstractFloat} = Float32, kwargs...)
    if (rg.orbit ÷ 360) * 360 == rg.orbit
        return ones(T, 1, 1, 1) # (1,1,1) for type stability
    end
    weight = parker_weight_fan_short(rg ; kwargs...)
    weight = reshape(weight, dims(rg)[1], 1, dims(rg)[3]) # (ns,1,na)
    return weight
end

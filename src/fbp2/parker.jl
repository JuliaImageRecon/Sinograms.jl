# fbp2/parker.jl
# Parker weighting for FBP reconstruction

export parker_weight
# using Sinograms: SinoGeom, SinoPar, SinoFan, SinoMoj


"""
    parker_weight(sg::SinoGeom ; T = Float32)
Compute Parker weighting for non-360° orbits.
Returns `Vector{T}` of length `sg.na`.
"""
function parker_weight(sg::SinoPar ; T::DataType = Float32)::Vector{T}
    orbit = abs(sg.orbit)
    na = sg.na
    ad = @. abs(sg.ad - sg.orbit_start)

    wt = ones(T, na)

    if (sg.orbit ÷ 180) * 180 == orbit
        return wt
    end

    if orbit < 180
        @warn("orbit $orbit < 180")
        return wt
    end

    orbit > 360 && throw("only 180 ≤ |orbit| ≤ 360 supported for Parker weighting")
    extra = orbit - 180 # extra beyond 180

    ii = ad .< extra
    wt[ii] = @. abs2(sin(ad[ii] / extra * π/2))
    ii = ad .≥ 180
    wt[ii] = @. abs2(sin((orbit - ad[ii]) / extra * π/2))
    wt *= orbit / 180 # trick because of the back-projector normalization
    return wt
end


function parker_weight_fan(na::Int, orbit::RealU; T::DataType = Float32)
    wt = ones(T, na)
    if (orbit ÷ 360) * 360 == orbit
        return wt
    end
    throw("todo: short-scan fan-beam Parker weighting not done")
    wt = zeros(T, na)
    return wt
end

parker_weight(sg::SinoFan; T::DataType = Float32)::Vector{T} =
    parker_weight_fan(sg.na, sg.orbit; T)


function parker_weight(sg::SinoMoj ; T::DataType = Float32)::Vector{T}
    orbit = abs(sg.orbit)
    na = sg.na
    ((sg.orbit ÷ 180) * 180 == orbit) ||
        throw("No Parker weighting for Mojette geometry with orbit=$orbit")
    wt = ones(T, na)
    return wt
end

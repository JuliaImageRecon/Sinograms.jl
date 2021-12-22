export fbp_ramp
export ramp_flat
export ramp_arc


"""
    h, nn = fbp_ramp(how, n, ds, dsd)
'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
This sampled band-limited approach avoids the aliasing that would be
caused by sampling the ramp directly in the frequency domain.

in
- `how::Symbol`                     `:arc` (3rd generation CT) or `:flat` (for parallel too)
- `n::Int`                          # of samples (must be even)
- `ds::RealU`                       sample spacing (in distance units, e.g., cm)
- `dsd::RealU`                      source-to-detector distance, for `:arc` case

out
- `h::AbstractMatrix{<:Number}`     samples of band-limited ramp filter

"""
fbp_ramp(::SinoPar, n::Int, ds::RealU, dsd::RealU) = ramp_flat(n, ds)
fbp_ramp(::SinoFanFlat, n::Int, ds::RealU, dsd::RealU) = ramp_flat(n, ds)
fbp_ramp(::SinoFanArc, n::Int, ds::RealU, dsd::RealU) = ramp_arc(n, ds, dsd)

function ramp_arc(n::Int, ds::RealU, dsd::RealU)

    isodd(n) && throw("n must be even")

    if n/2 * ds / dsd > 0.9 * pi/2
        println("angle is $(rad2deg(n/2 * ds / dsd)) degrees: too large")
        throw("physically impossible arc geometry")
    end

    nn = -(n÷2):(n÷2-1)
    h = zeros(size(nn))
    h[nn.==0] .= 1 / (4 * abs2(ds))
    odd = isodd.(nn)
    h[odd] = @.( -1 / abs2((pi * dsd * sin(nn[odd] * ds / dsd))))

    return h, nn 
end

function ramp_flat(n::Int, ds::RealU)

    isodd(n) && throw("n must be even")

    nn = -(n÷2):(n÷2-1)
    h = zeros(size(nn))
    h[n÷2+1] = 1 / 4
    odd = isodd.(nn)
    h[odd] .= -1 ./ abs2.(pi * nn[odd])
    h ./= abs2(ds)

    return h, nn
end


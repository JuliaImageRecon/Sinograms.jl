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
function fbp_ramp(how::Symbol, n::Int, ds::RealU, dsd::RealU)
    
    isodd(n) && throw("n must be even")

    # todo: replace if statements with multiple dispatch
    if how==:arc
        
	    h, nn = ramp_arc(n, ds, dsd)

    elseif how==:flat   
	    h, nn = ramp_flat(n, ds)

    else
        throw("bad fan type")
    end
    return h, nn
end


function ramp_arc(n::Int, ds::RealU, dsd::RealU)

    if n/2 * ds / dsd > 0.9 * pi/2
        println("angle is $(rad2deg(n/2 * ds / dsd)) degrees: too large")
        throw("physically impossible arc geometry")
    end

    nn = -(n÷2):(n÷2-1)
    h = zeros(size(nn))
    h[nn.==0] .= 1 / (4 * abs2(ds))
    odd = isodd.(nn)
    h[odd] .= abs2.(@.(-1 / (pi * dsd * sin(nn[odd] * ds / dsd))))

    return h, nn 
end


function ramp_flat(n::Int, ds::RealU)
    nn = -(n÷2):(n÷2-1)
    h = zeros(size(nn))
    h[n÷2+1] = 1 / 4
    odd = isodd.(nn)
    h[odd] .= -1 ./ abs2.(pi * nn[odd])
    h ./= abs2(ds)

    return h, nn
end


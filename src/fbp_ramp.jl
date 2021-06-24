export fbp_ramp

"""
    h, nn = fbp_ramp(how, n, ds, dsd)
'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
This sampled band-limited approach avoids the aliasing that would be
caused by sampling the ramp directly in the frequency domain.

in
- `how::Symbol`                     `:arc` (3rd generation CT) or `:flat` (for parallel too)
- `n::Int`                          # of samples (must be even)
- `ds::Real`                        sample spacing (in distance units, e.g., cm)
- `dsd::Real`                       source-to-detector distance, for `:arc` case

out
- `h::AbstractMatrix{<:Number}`     samples of band-limited ramp filter

"""
function fbp_ramp(how::Symbol, n::Int, ds::Real, dsd::Real)
    
    isodd(n) && throw("n must be even")

    if how==:arc
        
        #if nargin ~= 4, help(mfilename), error 'need 4 args', end
	    h, nn = ramp_arc(n, ds, dsd)
        

    elseif how==:flat
        
        
        #if nargin ~= 3 && ~isempty(dsd), warn('only need 3 args for flat'), end
	    h, nn = ramp_flat(n, ds)
        

    else
        throw("bad fan type")
    end
    return h, nn

end


function ramp_arc(n, ds, dsd)

    if n/2 * ds / dsd > 0.9 * pi/2
        println("angle is $(rad2deg(n/2 * ds / dsd)) degrees: too large")
        throw("physically impossible arc geometry")
    end

    nn = -(n/2):(n/2-1)
    h = zeros(size(nn))
    h[nn.==0] .= 1 / (4 * ds^2)
    odd = nn .% 2 .== 1
    h[odd] .= -1 ./ (pi .* dsd .* sin.(nn[odd] .* ds / dsd)).^2

    return h, nn 

end


function ramp_flat(n, ds)
    nn = -(n/2):(n/2-1)
    h = zeros(size(nn))
    h[n/2+1] .= 1 / 4
    odd = nn .% 2 .== 1
    h[odd] .= -1 ./ (pi .* nn[odd]).^2
    h = h ./ ds^2

    return h, nn

end


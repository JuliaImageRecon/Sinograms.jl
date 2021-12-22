export fbp2_sino_filter

using FFTW

"""
    sino, Hk, hh, nn = fbp2_sino_filter(how, sino; 
    ds=1, dsd=0, extra=0, npad=0, decon1=1, window=:none)

Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
Both parallel-beam and fan-beam tomographic geometries are supported.
This approach of sampling the band-limited ramp avoids the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

in
- `how::Symbol`                        `:arc` (3rd generation CT) or `:flat` (for parallel too)
- `sino::AbstractArray{<:Number}`        `[nb (L)]` sinograms

options
- `ds::RealU`                           sample spacing (in distance units, e.g., cm) (default: 1)
- `dsd::RealU`                          source-to-detector distance, for `:arc` case only. (default: Inf)
- `extra::Int`                          # of extra sinogram radial samples to keep (default: 0)
- `npad::Int`                           # of padded samples. (default: 0, means next power of 2)
- `decon1::Int`                         deconvolve effect of linear interpolator? (default: 1)
- `window::Symbol`                      default: `:none`

out
- `sino::AbstractArray{<:Number}`       filtered sinogram rows
- `Hk::AbstractVector{<:Number}`        apodized ramp filter frequency response
- `hn::AbstractVector{<:Number}`        samples of band-limited ramp filter
- `nn::AbstractVector{<:Number}`        [-np/2,...,np/2-1] vector for convenience

"""
function fbp2_sino_filter(
    sg::SinoGeom, 
    sino::AbstractMatrix{<:Number}; 
    ds::RealU=1, 
    dsd::RealU=Inf, 
    extra::Int=0, 
    npad::Int=0, 
    decon1::Int=1, 
    window::Symbol=:none
    )

    nb, na = size(sino)
    if npad==0
        npad = 2^ceil(Int64, log2(2*nb-1)) # padded size
    end

    sino = [sino; zeros(npad-nb,na)] # padded sinogram
    hn, nn = fbp_ramp(sg, npad, ds, dsd)

    reale = (x) -> (@assert x ≈ real(x); real(x))
    Hk = fft(fftshift(hn))
    Hk = reale(Hk)

    Hk = Hk .* fbp2_window(npad, window)

    Hk = ds * Hk # differential for discrete-space convolution vs integral


    #= linear interpolation is like blur with a triangular response,
    so we can compensate for this approximately in frequency domain =#
    if decon1 != 0
        Hk = Hk ./ fftshift(sinc.(nn / npad).^2)
    end
    
    sino = reale(ifft(fft(sino, 1) .* repeat(Hk, 1, na), 1)) # apply filter
    #NOTE: was fft(sino, [], 1) and ifft(..., [], 1) in matlab  
    #NOTE: still assertion error for reale 
    # todo: definitely use broadcast or map here.  should be no need to repeat
    # trick: possibly keep extra column(s) for zeros! 
    sino = sino[1:(nb+extra),:]
    sino[(nb+1):(nb+extra),:] .= 0
    sino = reshape(sino, nb, na)

    return sino, Hk, hn, nn
end
 
function fbp2_sino_filter(how::Symbol, sino::AbstractArray{<:Number}; kwargs...)
    return mapslices(sino -> fbp2_sino_filter(how, sino; kwargs...), sino, [1,2])
end


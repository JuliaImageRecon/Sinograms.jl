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
- `how::Symbol`                         `:arc` (3rd generation CT) or `:flat` (for parallel too)
- `sino::AbstractVector{<:Real}`        `[nb (L)]` sinograms

options
- `ds::Real=1`                          sample spacing (in distance units, e.g., cm) (default: 1)
- `dsd::Real=0`                         source-to-detector distance, for `:arc` case only.
- `extra::Int=0`                        # of extra sinogram radial samples to keep (default: 0)
- `npad::Int=0`                         # of padded samples. (default: 0, means next power of 2)
- `decon1::Int=1`                       deconvolve effect of linear interpolator? (default: 1)
- `window::Symbol=:none`

out
- `sino::AbstractMatrix{<:Number}`      filtered sinogram rows
- `Hk::AbstractMatrix{<:Number}`        apodized ramp filter frequency response
- `hn::AbstractMatrix{<:Number}`        samples of band-limited ramp filter
- `nn::AbstractMatrix{<:Number}`        [-np/2,...,np/2-1] vector for convenience

"""
function fbp2_sino_filter(how::Symbol, sino::AbstractVector{<:Real}; 
    ds::Real=1, dsd::Real, extra::Int=0, npad::Int=0, decon1::Int=1, window::Symbol=:none)


    dims = size(sino)
    sino = reshape(sino, dims[1], :)

    
    nb, na = size(sino)
    if npad==0
        npad = 2^ceil(log2(2*nb-1)) # padded size
    end

    sino = [sino; zeros(npad-nb,na)] # padded sinogram

    hn, nn = fbp_ramp(how, npad, ds, dsd)

    reale = (x) -> (@assert x ≈ real(x); real(x))

    Hk = reale.(fft(fftshift(hn)))

    Hk = Hk .* fbp2_window(npad, window)

    Hk = ds .* Hk # differential for discrete-space convolution vs integral

    

    #= linear interpolation is like blur with a triangular response,
    so we can compensate for this approximately in frequency domain =#
    if decon1 != 0
        Hk = Hk ./ fftshift(nufft_sinc(nn / npad).^2)
    end

    sino = ifft_sym( fft(sino, [], 1) .* repeat(Hk, [1 na]), [], 1) # apply filter

    # trick: possibly keep extra column(s) for zeros!
    sino = sino[1:(nb+extra),:]
    sino[(nb+1):(nb+extra),:] .= 0

    
    sino = reshape(sino, (size(sino, 1) dims[2:end]))

    return sino, Hk, hn, nn

end




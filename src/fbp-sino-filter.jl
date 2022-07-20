# fbp-sino-filter.jl

export fbp_sino_filter

using FFTW
#using Sinograms: SinoGeom, fbp_ramp, fbp2_window


"""
    sino, Hk, hh, nn = fbp_sino_filter(sg::SinoGeom, sino ;
        extra=0, npad=0, decon1=1, window=:none)

Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
Both parallel-beam and fan-beam tomographic geometries are supported.
This approach of sampling the band-limited ramp avoids the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

in
- `sg::SinoGeom`
- `sino::AbstractArray{<:Number}` : `[nb (L)]` sinograms

options
- `extra::Int` : # of extra sinogram radial samples to keep (default: 0)
- `npad::Int` : # of padded samples. (default: 0, means next power of 2)
- `decon1::Int` : deconvolve effect of linear interpolator? (default: 1)
- `window::Symbol` : default: `:none`

out
- `sino::AbstractArray`   filtered sinogram rows
- `Hk::Vector{<:Number}`  apodized ramp filter frequency response
- `hn::Vector{<:Number}`  samples of band-limited ramp filter
- `nn::AbstractVector`    `[-np/2,...,np/2-1]` vector for convenience

"""
function fbp_sino_filter(
    sg::SinoGeom,
    sino::AbstractArray{<:Number};
    extra::Int = 0,
    npad::Int = nextpow(2, sg.nb + 1),
    decon1::Int = 1,
    window::Symbol = :none,
)
@show sg.nb npad

    ds = sg.d

    sg.nb == size(sino,1) || throw("sinogram nb mismatch")
    sg.na == size(sino,2) || throw("sinogram na mismatch")
    nb = sb.nb
    na = sg.na
#   if npad == 0
#       npad = 2^ceil(Int64, log2(2*nb-1)) # padded size
#   end

    nb + extra > npad && throw("nb=$nb + extra=$extra > npad=$npad")

    dim = size(sino)
    dim[2] = npad
    sino = cat(dims=2, sino, zeros(dim .- size(sino)))
#   sino = [sino; zeros(npad-nb,na)] # padded sinogram
    hn, nn = fbp_ramp(sg, npad)

    reale = (x) -> (@assert x ≈ real(x); real(x))
    Hk = fft(fftshift(hn))
    Hk = reale(Hk)

    Hk .*= fbp2_window(npad, window)
    Hk = ds * Hk # differential for discrete-space convolution vs integral

    #= linear interpolation is like blur with a triangular response,
    so we can compensate for this approximately in frequency domain =#
    if decon1 != 0
        Hk ./= fftshift(sinc.(nn / npad).^2)
    end

    sino = ifft(fft(sino, 1) .* Hk, 1) # apply filter to each sinogram row
    sino = reale(sino) # todo: only if real data
    #NOTE: still assertion error for reale
    # todo: definitely use broadcast or map here.  should be no need to repeat
    # trick: possibly keep extra column(s) for zeros!
    sino = reshape(sino, npad, :)
    sino = sino[1:(nb+extra),:]
    sino[(nb+1):(nb+extra),:] .= zero(eltype(sino))
    sino = reshape(sino, nb, na, dim[3:end]...)

    return sino, Hk, hn, nn
end

#function fbp2_sino_filter(how::Symbol, sino::AbstractArray{<:Number}; kwargs...)
#    return mapslices(sino -> fbp2_sino_filter(how, sino; kwargs...), sino, [1,2])
#end

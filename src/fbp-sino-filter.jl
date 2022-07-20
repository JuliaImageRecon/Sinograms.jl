# fbp-sino-filter.jl

export fbp_sino_filter

using FFTW
#using Sinograms: SinoGeom, fbp_ramp, fbp2_window


"""
    sino, Hk, hh, nn = fbp_sino_filter(sg::SinoGeom, sino ;
        extra=0, npad=0, decon1=1, window=:none)

Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
Both parallel-beam and fan-beam tomographic geometries are supported.
This code samples the band-limited ramp to avoid the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

in
- `sg::SinoGeom`
- `sino::AbstractArray{<:Number}` : `[nb (L)]` sinograms

options
- `extra::Int` # of extra sinogram radial samples to keep (default: 0)
- `npad::Int` # of padded samples. (default: next power of 2)
- `decon1::Bool` deconvolve effect of linear interpolator? (default: true)
- `window::Symbol` apodizer; default: `:none`

out
- `sino::AbstractArray` sinogram with filtered rows
- `Hk::Vector` apodized ramp filter frequency response
- `hn::Vector` samples of band-limited ramp filter
- `nn::AbstractVector` `-(npad/2):(npad/2-1)` vector for convenience

"""
function fbp_sino_filter(
    sg::SinoGeom,
    sino::AbstractArray{<:Number};
    extra::Int = 0,
    npad::Int = nextpow(2, sg.nb + 1),
    decon1::Bool = true,
    window::Symbol = :none,
)
@show sg.nb npad

    ds = sg.d

    nb = sg.nb
    na = sg.na
    nb == size(sino,1) || throw("sinogram nb mismatch")
    na == size(sino,2) || throw("sinogram na mismatch")
#   if npad == 0
#       npad = 2^ceil(Int64, log2(2*nb-1)) # padded size
#   end
    nb + extra > npad && throw("nb=$nb + extra=$extra > npad=$npad")

    dimpadding = collect(size(sino))
    dimpadding[1] = npad - dimpadding[1]
#   tmp = dim - collect(size(sino))
@show dimpadding
    tmp = zeros(eltype(sino), dimpadding...)
    sino = cat(dims=1, sino, tmp)
@show size(sino)
#   sino = [sino; zeros(npad-nb,na)] # padded sinogram
    hn, nn = fbp_ramp(sg, npad)

    reale = (x) -> (@assert x ≈ real(x); real(x))
    Hk = fft(fftshift(hn))
    Hk = reale(Hk)

    Hk .*= fbp2_window(npad, window)
    Hk = ds * Hk # differential for discrete-space convolution vs integral

    # Linear interpolation is like blur with a triangular response,
    # so we can compensate for this approximately in frequency domain.
    if decon1 != 0
        Hk ./= fftshift(sinc.(nn / npad).^2)
    end

    is_real = isreal(sino)
    sino = ifft(fft(sino, 1) .* Hk, 1) # apply filter to each sinogram row
    if is_real
        sino = reale(sino)
    end

    # trick: possibly keep extra column(s) for zeros!
    sino = reshape(sino, npad, :) # (npad,…)
    sino = sino[1:(nb+extra),:]
    sino[(nb+1):(nb+extra),:] .= zero(eltype(sino))
    sino = reshape(sino, nb, na, dimpadding[3:end]...) # for >2D sinogram array

    return sino, Hk, hn, nn
end

#=
function fbp_sino_filter(how::Symbol, sino::AbstractArray{<:Number}; kwargs...)
    return mapslices(sino -> fbp2_sino_filter(how, sino; kwargs...), sino, [1,2])
end
=#

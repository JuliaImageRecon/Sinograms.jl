# fbp2/sino-filter.jl

export fbp_filter, fbp_sino_filter

using FFTW
#using Sinograms: SinoGeom, SinoPar, Window, fbp_ramp, fbp_window


function _reale(x)
    (x ≈ real(x)) || @warn("x not real $(maximum(abs, x-real(x))/maximum(abs,x))")
    return real(x)
end


"""
    Hk = fbp_filter(sg::SinoGeom ;
        npad=0, ds::RealU = sg.d, decon1::Bool=true, window=Window())

Compute frequency response of ramp-like filter
used for 2D FBP image reconstruction.
Both parallel-beam and fan-beam tomographic geometries are supported.
This code samples the band-limited ramp to avoid the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

in
- `sg::SinoGeom`

options
- `npad::Int` # of padded samples. (default: next power of 2)
- `decon1::Bool` deconvolve effect of linear interpolator? (default: true)
- `window::Window` apodizer; default: `Window()`

out
- `Hk::Vector` apodized ramp filter frequency response
"""
function fbp_filter(
    sg::SinoGeom{Td} = SinoPar() ;
    npad::Int = nextpow(2, sg.nb + 1),
    ds::Td = sg.d,
    decon1::Bool = true,
    window::Window = Window(),
) where {Td <: RealU}

#   U = eltype(1 / oneunit(Td))
    hn, nn = fbp_ramp(sg, npad)

    unit = oneunit(eltype(hn)) # handle units
    Hk = unit * fft(fftshift(hn / unit))
    Hk = _reale(Hk)

    Hk .*= fbp_window(window, npad)
    Hk = ds * Hk # differential for discrete-space convolution vs integral

    # Linear interpolation is like blur with a triangular response,
    # so we can compensate for this approximately in frequency domain.
    if decon1
        Hk ./= fftshift(sinc.(nn / npad).^2)
    end

    return Hk#::Vector{U}
end


"""
    sino, Hk, hh, nn = fbp_sino_filter(sg::SinoGeom, sino ;
        extra=0, npad=0, decon1=1, window=Window())

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
- `window::Window` apodizer; default: `Window()`

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
    window::Window = Window(),
)

    ds = sg.d
    nb = sg.nb
    na = sg.na
    nb == size(sino,1) || throw("sinogram nb mismatch")
    na == size(sino,2) || throw("sinogram na mismatch")
    nb + extra > npad && throw("nb=$nb + extra=$extra > npad=$npad")

    dimpadding = collect(size(sino))
    dimpadding[1] = npad - dimpadding[1]
    tmp = zeros(eltype(sino), dimpadding...)
    sino = cat(dims=1, sino, tmp) # padded sinogram
    hn, nn = fbp_ramp(sg, npad)

    unit = oneunit(eltype(hn)) # handle units
    Hk = unit * fft(fftshift(hn / unit))
    Hk = _reale(Hk)

    Hk .*= fbp_window(window, npad)
    Hk = ds * Hk # differential for discrete-space convolution vs integral

    # Linear interpolation is like blur with a triangular response,
    # so we can compensate for this approximately in frequency domain.
    if decon1 != 0
        Hk ./= fftshift(sinc.(nn / npad).^2)
    end

    is_real = isreal(sino)
    unit = oneunit(eltype(sino)) # handle units
    sino = unit * ifft(fft(sino/unit, 1) .* Hk, 1) # apply filter to each sinogram row
    if is_real
        sino = _reale(sino)
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
    return mapslices(sino -> fbp_sino_filter(how, sino; kwargs...), sino, [1,2])
end
=#

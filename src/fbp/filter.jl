# fbp/filter.jl

export fbp_filter, fbp_sino_filter

using FFTW: fft, ifft, fftshift

#using Sinograms: SinoGeom, SinoPar, Window, fbp_ramp, fbp_window, _reale


"""
    Hk = fbp_filter(rg::RayGeom ;
        npad=0, ds::RealU = rg.d, decon1::Bool=true, window=Window())

Compute frequency response of ramp-like filter
used for FBP image reconstruction.
Supports parallel-beam and fan-beam tomographic geometries in 2D and 3D.
This code samples the band-limited ramp to avoid the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

# in
- `rg::RayGeom`

# option
- `npad::Int` # of padded samples. (default: next power of 2)
- `ds::Td` detector sample spacing (default from `st`)
- `decon1::Bool` deconvolve effect of linear interpolator? (default: true)
- `window::Window` apodizer; default: `Window()`

# out
- `Hk::Vector` apodized ramp filter frequency response
"""
function fbp_filter(
    rg::RayGeom{Td} = SinoPar(),
    ;
    npad::Int = nextpow(2, rg.nb + 1),
    ds::Td = _ds(rg),
    decon1::Bool = true,
    window::Window = Window(),
) where {Td <: RealU}

    hn, nn = fbp_ramp(rg, npad)

    unit = oneunit(eltype(hn)) # handle units
    Hk = unit * fft(fftshift(hn / unit))
    Hk = _reale(Hk)

    Hk .*= fbp_window(window, npad)
    Hk = ds * Hk # differential for discrete-space convolution vs integral

    # Linear interpolation is like blur with a triangular response,
    # so we can approximately compensate for this effect in the frequency domain.
    if decon1
        Hk ./= fftshift(sinc.(nn / npad).^2)
    end

    return Hk
end


"""
    fft_filter(data::Array, filter::Vector [, dim])

Apply filter to selected dimensions of array `data` using FFT.

# in
- `data::AbstractArray{<:Number} (n, (L))`
- `filter::AbstractVector (n)` apodized ramp filter frequency response

# option
- `dim` non-singleton dimensions of `filter` (typically `1`)

# out
- `out::AbstractArray` data filtered along dimension `dim`

If the input data is real,
so will be the output;
this assumes `filter` frequency response
has appropriate symmetry.
"""
fft_filter

function _fft_filter(
    data::AbstractArray{<:Number},
    filter::AbstractArray,
    dim = findall(!=(1), size(filter)),
)
    return ifft(fft(data, dim) .* filter, dim)
end

function fft_filter(
    data::AbstractArray{<:Complex},
    args...
)
    return _fft_filter(data, args...)
end

function fft_filter(
    data::AbstractArray{<:Real},
    args...
)
    return _reale(_fft_filter(data, args...))
end


"""
    sino = fbp_sino_filter(sino::Array, filter::Vector ; extra=0)

Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
Supports both parallel-beam and fan-beam tomographic geometries in 2D and 3D.

# in
- `sino::AbstractArray{<:Number}` `[nb (L)]` sinograms
- `filter::AbstractVector` `(npad ≥ nb)` apodized ramp filter frequency response

# option
- `extra::Int` # of extra sinogram radial samples to keep (default: 0)
- `npad::Int` # of padded samples. (default: next power of 2)

# out
- `sino::AbstractArray` sinogram with filtered rows
"""
function fbp_sino_filter(
    sino::AbstractArray{Ts},
    filter::AbstractVector{<:Number},
    ;
    extra::Int = 0,
) where {Ts <: Number}

    npad = length(filter)
    nb = size(sino, 1)
    nb + extra > npad && error("nb=$nb + extra=$extra > npad=$npad")

    dimpadded = (npad, size(sino)[2:end]...)
    tmp = zeros(Ts, dimpadded)
    selectdim(tmp, 1, 1:size(sino,1)) .= sino
    sino = tmp # padded sinogram

    sino = fft_filter(sino, filter) # apply filter to each sinogram row

    # trick: possibly keep extra column(s) for zeros!
    sino = reshape(sino, npad, :) # (npad,…)
    sino = @view sino[1:(nb+extra), :]
    sino[(nb+1):(nb+extra), :] .= zero(eltype(sino))
    sino = reshape(sino, nb, :, dimpadded[3:end]...) # for >2D sinogram array

    return sino
end

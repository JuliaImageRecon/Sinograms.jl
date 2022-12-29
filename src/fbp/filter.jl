# fbp/filter.jl

export fbp_filter, fbp_sino_filter

using FFTW
#using Sinograms: SinoGeom, SinoPar, Window, fbp_ramp, fbp_window, _reale


"""
    Hk = fbp_filter(rg::RayGeom ;
        npad=0, ds::RealU = rg.d, decon1::Bool=true, window=Window())

Compute frequency response of ramp-like filter
used for FBP image reconstruction.
Supports parallel-beam and fan-beam tomographic geometries in 2D and 3D.
This code samples the band-limited ramp to avoid the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

in
- `rg::RayGeom`

options
- `npad::Int` # of padded samples. (default: next power of 2)
- `ds::Td` detector sample spacing (default from `st`)
- `decon1::Bool` deconvolve effect of linear interpolator? (default: true)
- `window::Window` apodizer; default: `Window()`

out
- `Hk::Vector` apodized ramp filter frequency response
"""
function fbp_filter(
    rg::RayGeom{Td} = SinoPar() ;
    npad::Int = nextpow(2, rg.nb + 1),
    ds::Td = _ds(rg),
    decon1::Bool = true,
    window::Window = Window(),
) where {Td <: RealU}

#   U = eltype(1 / oneunit(Td))
    hn, nn = fbp_ramp(rg, npad)

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
    sino = fbp_sino_filter(sino::Array, filter::Vector ; extra=0)

Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
Supports both parallel-beam and fan-beam tomographic geometries in 2D and 3D.

in
- `sino::AbstractArray{<:Number}` `[nb (L)]` sinograms
- `filter::AbstractVector` `(npad ≥ nb)` apodized ramp filter frequency response

options
- `extra::Int` # of extra sinogram radial samples to keep (default: 0)
- `npad::Int` # of padded samples. (default: next power of 2)

out
- `sino::AbstractArray` sinogram with filtered rows
"""
function fbp_sino_filter(
    sino::AbstractArray{Ts, N},
    filter::AbstractVector{Tf} ;
    extra::Int = 0,
) where {N, Ts <: Number, Tf <: Number}

    npad = length(filter)
    nb = size(sino,1)
    nb + extra > npad && throw("nb=$nb + extra=$extra > npad=$npad")

    dimpadding = collect(size(sino))
    dimpadding[1] = npad - dimpadding[1]
    tmp = zeros(Ts, dimpadding...)
    sino = cat(dims=1, sino, tmp) # padded sinogram

    is_real = isreal(sino)

    # handle fft with units
    unit_s = oneunit(Ts)
    unit_f = oneunit(Tf)
    # apply filter to each sinogram row
    sino = (unit_s * unit_f) * ifft(fft(sino/unit_s, 1) .* (filter/unit_f), 1)

    if is_real
        sino = _reale(sino)
    end

    # trick: possibly keep extra column(s) for zeros!
    sino = reshape(sino, npad, :) # (npad,…)
    sino = sino[1:(nb+extra), :]
    sino[(nb+1):(nb+extra), :] .= zero(eltype(sino))
    sino = reshape(sino, nb, :, dimpadding[3:end]...) # for >2D sinogram array

    U = eltype(unit_s * unit_f)
    return sino::Array{U, N}
end


#=
function fbp_sino_filter(
    aa::AbstractArray{Ts,N},
    filter::AbstractVector{Tf} ;
    kwargs...,
) where {N, Ts <: Number, Tf <: Number}
    fun(sino) = fbp_sino_filter(sino, filter; kwargs...)
    unit_s = oneunit(Ts)
    unit_f = oneunit(Tf)
    U = eltype(unit_s * unit_f)
    out = mapslices(fun, aa, dims=[1,2])
    return out::Array{U,N}
end
=#


#= old version with window
"""
    sino, Hk, hh, nn = fbp_sino_filter(rg::SinoGeom, sino ;
        extra=0, npad=0, decon1=1, window=Window())

Apply ramp-like filters to sinogram(s) for 2D FBP image reconstruction.
Both parallel-beam and fan-beam tomographic geometries are supported.
This code samples the band-limited ramp to avoid the aliasing that
would be caused by sampling the ramp directly in the frequency domain.

in
- `rg::SinoGeom`
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
    rg::SinoGeom,
    sino::AbstractArray{<:Number};
    extra::Int = 0,
    npad::Int = nextpow(2, rg.nb + 1),
    decon1::Bool = true,
    window::Window = Window(),
)

    ds = rg.d
    nb = rg.nb
    na = rg.na
    nb == size(sino,1) || throw("sinogram nb mismatch")
    na == size(sino,2) || throw("sinogram na mismatch")
    nb + extra > npad && throw("nb=$nb + extra=$extra > npad=$npad")

    dimpadding = collect(size(sino))
    dimpadding[1] = npad - dimpadding[1]
    tmp = zeros(eltype(sino), dimpadding...)
    sino = cat(dims=1, sino, tmp) # padded sinogram
    hn, nn = fbp_ramp(rg, npad)

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
=#

#=
function fbp_sino_filter(how::Symbol, sino::AbstractArray{<:Number}; kwargs...)
    return mapslices(sino -> fbp_sino_filter(how, sino; kwargs...), sino, [1,2])
end
=#

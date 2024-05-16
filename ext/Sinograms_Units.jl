# Sinograms_Units.jl
# Support data with units iff user has loaded Unitful

module Sinograms_Units

import Unitful: °, rad, NoDims, unit
import Unitful: Units, Quantity, convfact, ustrip
import Unitful: uconvert
import Unitful
import FFTW #: fft, ifft

#=
https://github.com/PainterQubits/Unitful.jl/issues/375
=#
function uconvert(a::Units, x::Quantity{T,D,U}) where {T<:AbstractFloat,D,U}
    return Quantity(x.val * T(convfact(a, U())), a)
end


"""
    to_radians(angles::AbstractArray{Unitful.Quantity})

Convert `Unitful` quantity array to radians.
"""
function to_radians(aa::AbstractArray{<: Unitful.Quantity{T}}) where {T <: AbstractFloat}
    U = eltype(aa)
    c = rad(oneunit(U)) / oneunit(U)
    return aa * c
end


_unit_precision(x::Unitful.Quantity{T}) where {T <: Number} = "Unit{$T}"


# generic unitless linear operation applied to data with units
function _linear_fun(fun::Function,
    x::AbstractArray{<: Unitful.Quantity},
    args...
)
    x0 = ustrip(x) # unitless view into x data
    y0 = fun(x0, args...)
    u = unit(eltype(x))
    Tu = typeof(one(eltype(y0)) * unit(eltype(x)))
    return reinterpret(Tu, y0)
end

# fft for data with units
FFTW.fft(x::AbstractArray{<: Unitful.Quantity}, args...) =
    _linear_fun(FFTW.fft, x, args...)
FFTW.ifft(x::AbstractArray{<: Unitful.Quantity}, args...) =
    _linear_fun(FFTW.ifft, x, args...)


fft_filter(data::AbstractArray{<: Unitful.Quantity{<:Complex}}, args...) =
    _fft_filter(data, args...)

fft_filter(data::AbstractArray{<: Unitful.Quantity{<:Real}}, args...) =
    _reale(_fft_filter(data, args...))

end # module

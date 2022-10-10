# units.jl

using .Unitful: °, rad, NoDims

#=
https://github.com/PainterQubits/Unitful.jl/issues/375
=#
using .Unitful: Units, Quantity, convfact
import .Unitful: uconvert
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

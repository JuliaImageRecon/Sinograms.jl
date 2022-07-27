# units.jl

using .Unitful: °, rad


"""
    to_radians(angle::Unitful.Quantity)
    to_radians(angles::AbstractArray{Unitful.Quantity})

Convert `Unitful` quantity to radians.
"""
#to_radians(angle::Unitful.Quantity) = convert(eltype(1.0rad), angle)
#to_radians(aa::AbstractArray{T}) where {T <: Unitful.Quantity} = convert(eltype(1.0rad), angle)
#to_radians(aa::AbstractArray{<: Unitful.Quantity{T}}) where {T <: Real} =
#    aa * ((one(Float32)*rad) / oneunit(eltype(aa)))
to_radians(aa::AbstractArray{<: Unitful.Quantity{T}}) where {T <: AbstractFloat} =
    aa * ((one(T)*rad) / oneunit(eltype(aa)))

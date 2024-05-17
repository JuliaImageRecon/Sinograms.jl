# unit.jl

export to_radians

"""
    to_radians(angles::AbstractArray{<:AbstractFloat})
When Unitful package not loaded,
assume `angles` are in degrees and convert to radians.
"""
to_radians(aa::AbstractArray{T}) where {T <: AbstractFloat} = aa * T(deg2rad(1))

_unit_precision(x::T) where {T <: Number} = T

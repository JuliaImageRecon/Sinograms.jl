# units.jl

using .Unitful: °, rad

# see Sinograms.jl
#_worktype(::Unitful.Quantity{T}) where {T <: AbstractFloat} = T
#_worktype(::Unitful.Quantity{T}) where {T <: Real} = Float32


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


#=
function sino_geom_taufun(sg::SinoGeom{<:Unitful.Quantity{T}}, x, y) where {T <: AbstractFloat}
    size(x) != size(y) && throw("bad x,y size")
    out = _sino_geom_taufun(sg, vec(x), vec(y))
@show T eltype(out)
    return out::Matrix{T}
end
=#

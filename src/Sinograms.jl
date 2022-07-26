"""
    Sinograms module
"""
module Sinograms

using Requires: @require

const RealU = Number # Real or Unitful

include("fbp/window.jl")

include("fbp2/sino-geom.jl")
include("fbp2/parker.jl")
include("fbp2/ramp.jl")
include("fbp2/sino-filter.jl")
include("fbp2/sino-weight.jl")
include("fbp2/back-par.jl")
include("fbp2/back-fan.jl")

include("fbp2.jl")
#include("fbp2_back.jl")
#include("fbp2_back_fan.jl")

include("fbp-par.jl")

include("sys2/zwart_powell.jl")


"""
    to_radians(angle::Real)
    to_radians(angles::AbstractArray{<:Real})
If no Unitful package loaded, assume `angle` is in degrees and convert to radians.
"""
#to_radians(angle::T) where {T <: AbstractFloat} = T(deg2rad(angle))
#to_radians(angle::Real) = deg2rad(Float32(angle))
to_radians(aa::AbstractArray{T}) where {T <: AbstractFloat} = aa * T(deg2rad(1))
#to_radians(aa::AbstractArray{T}) where {T <: Real} = aa * deg2rad(1f0) # Float32


# working precision should be at least Float32 (see also units.jl)
#_worktype(::T) where {T <: AbstractFloat} = T
#_worktype(::T) where {T <: Real} = Float32


# support Plots iff user has loaded that package
function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("fbp2/sino-plot.jl")
    @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" include("units.jl")
end


end # module

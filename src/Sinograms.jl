"""
    Sinograms module
"""
module Sinograms

using Requires: @require

const RealU = Number # Real or Unitful

include("fbp/reale.jl")
include("fbp/window.jl")

include("fbp2/sino-geom.jl")
include("fbp2/parker.jl")
include("fbp2/ramp.jl")
include("fbp2/sino-filter.jl")
include("fbp2/sino-weight.jl")
include("fbp2/back-par.jl")
include("fbp2/back-fan.jl")
include("fbp2/plan.jl")
include("fbp2/fbp.jl")

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
_unit_precision(x::T) where {T <: Number} = T


# support Plots (with units) iff user has loaded the relevant packages
function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
     @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" begin
      @require UnitfulRecipes = "42071c24-d89e-48dd-8a24-8a12d9b8861f" begin
       include("fbp2/sino-plot.jl")
      end
     end
    end

    @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" begin
        include("units.jl")
    end
end



end # module

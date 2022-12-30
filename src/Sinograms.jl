"""
    Sinograms module
"""
module Sinograms

using Requires: @require

const RealU = Number # Real or Unitful

include("geom/util.jl")
include("geom/ct-source.jl")
include("geom/type2.jl")
include("geom/type3.jl")
include("geom/sino-geom.jl")
include("geom/ct-geom.jl")
include("geom/common.jl")
include("geom/tau.jl")

include("fbp/reale.jl")
include("fbp/window.jl")

include("fbp/d-angle.jl")
include("fbp/ramp.jl")
include("fbp/filter.jl")
include("fbp/parker.jl")

include("fbp2/sino-weight.jl")
include("fbp2/back-par.jl")
include("fbp2/back-fan.jl")
include("fbp2/plan2.jl")
include("fbp2/fbp.jl")

include("fbp2/fbp-par.jl")
#include("fbp2/fbp-axis.jl")

include("fbp3/plan3.jl")
include("fbp3/cb_arc_to_par.jl")
include("fbp3/cb_flat_to_par.jl")
include("fbp3/cbct-back.jl")
include("fbp3/fdk.jl")

include("sys2/footprint.jl")
include("sys2/zwart_powell.jl")


"""
    to_radians(angles::AbstractArray{<:AbstractFloat})
When Unitful package not loaded,
assume `angles` are in degrees and convert to radians.
"""
to_radians(aa::AbstractArray{T}) where {T <: AbstractFloat} = aa * T(deg2rad(1))

_unit_precision(x::T) where {T <: Number} = T


# support Plots (with units) iff user has loaded the relevant packages
function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
       include("fbp2/sino-plot.jl")
       include("fbp3/ct-plot.jl")
    end

    @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" begin
        include("units.jl")
    end
end


end # module

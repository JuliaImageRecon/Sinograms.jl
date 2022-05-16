"""
    Sinograms module
"""
module Sinograms

const RealU = Number # Real or Unitful

include("fbp2.jl")
include("fbp2_sino_filter.jl")
include("fbp2_back.jl")
include("fbp2_back_fan.jl")
include("fbp_ramp.jl")
include("fbp2_window.jl")
include("fbp2_sino_weight.jl")

include("fbp-par.jl")

include("bdd_2d.jl")

end # module

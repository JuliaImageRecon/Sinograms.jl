"""
    Sinograms module
"""
module Sinograms

const RealU = Number # temporary

struct test
    num::Int
end


include("fbp2.jl")
include("fbp2_sino_filter.jl")
include("fbp2_back.jl")
include("fbp2_back_fan.jl")
include("fbp_ramp.jl")
include("fbp2_window.jl")



end # module

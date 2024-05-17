"""
    Sinograms module
"""
module Sinograms

const RealU = Number # Real or Unitful

include("unit.jl")
include("exts.jl")

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
include("sys2/bdd_2d.jl")

end # module

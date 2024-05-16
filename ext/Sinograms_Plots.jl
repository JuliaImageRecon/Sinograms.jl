# Sinograms_Plots.jl
# support plots (possibly with units) iff user has loaded Plots

module Sinograms_Plots

include("sino-plot.jl") # fbp2
include("ct-plot.jl") # fbp3

end # module

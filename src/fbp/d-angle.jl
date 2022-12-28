#=
fbp/d-angle.jl
=#

# scale projections by dβ for Riemann-like integration
function _view_weights(ar::AbstractVector{<:RealU})
    da = diff(ar) # typically 2π/na for 360° orbit
#   da1 = betas[begin] - (ar[end] - 2π)
    da = [da[1]; da] # todo
#   return reshape(da, 1, 1, :) # todo
    return da / 2
end

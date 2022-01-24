#=
sino-plot.jl
Show sinogram geometries.
2022-01-23, Jeff Fessler
=#

using Sinograms: RealU
using .Plots: scatter

"""
    scatter(sg::SinoGeom)

Make a scatter plot of the `(r, ϕ)` sample locations for all rays
"""
function scatter(sg::SinoGeom; kwargs...)
    r, phi = rays(sg)
    ylim = [min(0, rad2deg(minimum(phi))), max(360, rad2deg(maximum(phi)))]
    rmax = ceil(maximum(abs.(r))/10, digits=0)*10
    scatter(r, rad2deg.(phi), label="", markersize=1, markerstrokecolor=:auto,
        markershape=:circle, linewidth=0,
        ylim = ylim, xlim = [-1,1]*rmax, # ylabel="ϕ",
        title="$(typeof(sg))", xtick=(-1:1)*rmax, ytick=[0,360])
end


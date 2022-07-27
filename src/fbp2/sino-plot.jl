#=
sino-plot.jl
Show sinogram geometries.
2022-01-23, Jeff Fessler
=#

export sino_plot_rays

using Sinograms: RealU
using .Plots: scatter


"""
    sino_plot_rays(sg::SinoGeom ; kwargs...)

Make a scatter plot of the `(r, ϕ)` sample locations for all rays.
Requires `Plots`.
"""
function sino_plot_rays(sg::SinoGeom; kwargs...)
    r, phi = rays(sg)
    ad = rad2deg.(phi)
    unit_a = oneunit(eltype(ad))
    ylims = (min(zero(unit_a), minimum(ad)),
             max(360*unit_a, maximum(ad)))
    unit_r = oneunit(eltype(r))
    rmax = ceil(maximum(abs, r/unit_r)/10, digits=0) * 10 * unit_r
    xlims = (-1,1) .* rmax
    scatter(
        r, ad ;
        label="",
        markersize=1, markerstrokecolor=:auto,
        markershape=:circle, linewidth=0,
        ylims, xlims, ylabel="ϕ",
        xticks=(-1:1)*rmax, yticks=[0,360] * unit_a,
        title="$(typeof(sg))",
        kwargs...
    )
end

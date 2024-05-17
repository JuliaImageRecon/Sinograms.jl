#=
fbp2/sino-plot.jl
Show 2D sinogram geometries.
2022-01-23, Jeff Fessler
=#

import Sinograms: sino_plot_rays, sino_geom_plot! # extended here

using Sinograms: SinoGeom, RealU, rays, SinoPar, SinoFan, SinoMoj
using Sinograms: _ar, _rfov, _dso, _xds, _yds, _d_moj
using Plots: scatter, plot, plot!, default, xlims!
using ImageGeoms: ImageGeom


"""
    sino_plot_rays(rg::SinoGeom ; kwargs...)

Make a scatter plot of the `(r, ϕ)` sample locations for all rays.
Requires `Plots`.
"""
function sino_plot_rays(rg::SinoGeom; kwargs...)
#   r, phi = rays(rg)
    i = rays(rg)
    r = [i[1] for i in i]
    ϕ = [i[2] for i in i]
    ad = rad2deg.(ϕ)
    unit_a = oneunit(eltype(ad))
    ylims = (min(zero(unit_a), minimum(ad)),
             max(360*unit_a, maximum(ad)))
    unit_r = oneunit(eltype(r))
    rmax = ceil(maximum(abs, r/unit_r)/10, digits=0) * 10 * unit_r
    xlims = (-1,1) .* rmax
    scatter(
        r, ad ;
        markersize=1, markerstrokecolor=:auto,
        markershape=:circle, linewidth=0, label="",
        ylims, xlims, xlabel = "r", ylabel="ϕ",
        xticks = (-1:1)*rmax, yticks = [0,360] * unit_a,
        title = nameof(typeof(rg)),
        kwargs...
    )
end


# helpers

_round(x; kwargs...) = oneunit(x) * round(x / oneunit(x); kwargs...)

scat!(args... ; kwargs...) =
    plot!(args... ; kwargs..., markerstrokecolor=:auto, linewidth=0)


function sino_geom_plot_ig!(
    rfov::RealU ;
    ig::Union{Nothing,ImageGeom} = nothing,
)

    default(label="")

    xmax = rfov; xmin = -xmax; (ymin,ymax) = (xmin,xmax)
    if !isnothing(ig)
#       plot!(jim(ig.x, ig.y, ig.mask[:,:,1], clim=(0,1))) # todo: jim!
        x, y = axes(ig)
        xmin = minimum(x); xmax = maximum(x)
        ymin = minimum(y); ymax = maximum(y)
    end

    xticks = _round.([xmin, zero(xmin), xmax]; digits = 0)
    yticks = _round.([ymin, zero(ymin), ymax]; digits = 2)
    plot!(
        [xmax, xmin, xmin, xmax, xmax],
        [ymax, ymax, ymin, ymin, ymax] ;
        xlabel = "x", ylabel = "y",
        color = :green,
        xticks,
        yticks,
    )

    scat!([zero(xmin)], [zero(xmin)], marker=:circle)
    θ = range(0, 2π, 401)
    plot!(rfov * cos.(θ), rfov * sin.(θ), color=:magenta) # rfov circle
    plot!(aspect_ratio = 1)
end


function sino_geom_plot_fan!(
    ar,
    rfov,
    dso,
    xds,
    yds,
)

    x0 = zero(dso)
    y0 = dso
    t = range(0, 2π, 100)
    rot = ar[1]
    rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
    p0 = rot * [x0; y0]
    pd = rot * [xds'; yds'] # detector points

    xlims = (-1,1) .* dso
    ylims = (-1,1) .* dso
    tmp = ar .+ π/2 # trick: angle beta defined ccw from y axis
    scat!([p0[1]], [p0[2]], color=:blue, marker=:square, # source
        ; xlims, ylims, label="")
    plot!(dso * cos.(t), dso * sin.(t), color=:cyan, label="") # source circle
    scat!(dso * cos.(tmp), dso * sin.(tmp),
        color=:blue, marker=:circle, markersize=2, label="") # source points
    scat!(vec(pd[1,:]), vec(pd[2,:]), marker=:circle,
        color=:orange, markersize=1, label="") # detectors

    plot!([pd[1,1], p0[1], pd[1,end]], [pd[2,1], p0[2], pd[2,end]],
        color=:red, label="")
end


"""
    sino_geom_plot!(rg)
Picture of the source position / detector geometry.
"""
function sino_geom_plot!(rg::SinoGeom)

    plot!(; title = nameof(typeof(rg)))

#=
    if rg isa SinoPar
    end
=#

    if rg isa SinoFan
        sino_geom_plot_fan!(
            _ar(rg),
            _rfov(rg),
            _dso(rg),
            _xds(rg),
            _yds(rg),
        )

    elseif rg isa SinoMoj
        θ = range(0, 2π, 100)
        d_moj = _d_moj(rg)
        rphi = rg.nb/2 * d_moj.(θ)
        plot!(rphi .* cos.(θ), rphi .* sin.(θ), color=:blue, label="")
    #   rmax = maximum(_s(rg))
    #   axis([-1 1 -1 1] * max([rmax ig.fov/2]) * 1.1)
    end

    return plot!(aspect_ratio=1)
end


"""
    sino_geom_plot!(rg, ig)
Picture of the source position / detector geometry.
"""
function sino_geom_plot!(
    rg::SinoGeom,
    ig::ImageGeom;
)
    sino_geom_plot_ig!(_rfov(rg); ig)
    sino_geom_plot!(rg)
end

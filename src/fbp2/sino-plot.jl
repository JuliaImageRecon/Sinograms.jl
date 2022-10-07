#=
fbp2/sino-plot.jl
Show 2D sinogram geometries.
2022-01-23, Jeff Fessler
=#

export sino_plot_rays, sino_geom_plot!

#using Sinograms: SinoGeom, RealU
using .Plots: scatter, plot, plot!, default, xlims!


function _title(st::Union{SinoGeom,CtGeom})
    title = "$(typeof(st))"
    return title[1:findfirst('{', title)-1]
end


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
        markersize=1, markerstrokecolor=:auto,
        markershape=:circle, linewidth=0,
        ylims, xlims, xlabel = "r", ylabel="ϕ",
        xticks = (-1:1)*rmax, yticks = [0,360] * unit_a,
        title = _title(sg), label="",
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
    θ = LinRange(0, 2π, 401)
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
    t = LinRange(0, 2π, 100)
    rot = ar[1]
    rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
    p0 = rot * [x0; y0]
    pd = rot * [xds'; yds'] # detector points

    xlims = (-1,1) .* dso
    ylims = (-1,1) .* dso
    tmp = ar .+ π/2 # trick: angle beta defined ccw from y axis
    scat!([p0[1]], [p0[2]], color=:blue, marker=:square, # source
        ; xlims, ylims)
    plot!(dso * cos.(t), dso * sin.(t), color=:cyan) # source circle
    scat!(dso * cos.(tmp), dso * sin.(tmp),
        color=:blue, marker=:circle, markersize=2) # source points
    scat!(vec(pd[1,:]), vec(pd[2,:]), marker=:circle,
        color=:orange, markersize=1) # detectors

    plot!([pd[1,1], p0[1], pd[1,end]], [pd[2,1], p0[2], pd[2,end]],
        color=:red, label="")
end


"""
    sino_geom_plot!(sg)
Picture of the source position / detector geometry.
"""
function sino_geom_plot!(sg::SinoGeom)

    plot!(; title = _title(sg))

#=
    if sg isa SinoPar
    end
=#

    if sg isa SinoFan
        sino_geom_plot_fan!(
            sg.ar,
            sg.rfov,
            sg.dso,
            sg.xds,
            sg.yds,
        )

    elseif sg isa SinoMoj
        θ = LinRange(0, 2π, 100)
        rphi = sg.nb/2 * sg.d_moj.(θ)
        plot!(rphi .* cos.(θ), rphi .* sin.(θ), color=:blue, label="")
    #   rmax = maximum(sg.s)
    #   axis([-1 1 -1 1] * max([rmax ig.fov/2]) * 1.1)
    end

    return plot!(aspect_ratio=1)
end


"""
    sino_geom_plot!(sg, ig)
Picture of the source position / detector geometry.
"""
function sino_geom_plot!(
    sg::SinoGeom,
    ig::ImageGeom;
)
    sino_geom_plot_ig!(sg.rfov; ig)
    sino_geom_plot!(sg)
end

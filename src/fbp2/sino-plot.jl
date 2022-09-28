#=
sino-plot.jl
Show sinogram geometries.
2022-01-23, Jeff Fessler
=#

export sino_plot_rays, sino_geom_plot!

#using Sinograms: SinoGeom, RealU
using .Plots: scatter, plot!, default
using .Unitful
using .UnitfulRecipes


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
        ylims, xlims, ylabel="ϕ",
        xticks=(-1:1)*rmax, yticks=[0,360] * unit_a,
        title="$(typeof(sg))", label="",
        kwargs...
    )
end


"""
    sino_geom_plot!(sg ; ig)
Picture of the source position / detector geometry.
"""
function sino_geom_plot!(
    sg::SinoGeom ;
    ig::Union{Nothing,ImageGeom} = nothing,
)
    u = oneunit(sg.rfov) # hack to avoid unit issues

    default(label="")
    scat!(args... ; kwargs...) =
        plot!(args... ; kwargs..., markerstrokecolor=:auto, linewidth=0)

    xmax = sg.rfov / u; xmin = -xmax; (ymin,ymax) = (xmin,xmax)
    if !isnothing(ig)
#       plot!(jim(ig.x, ig.y, ig.mask[:,:,1], clim=(0,1))) # todo: jim!
        x, y = axes(ig)
        x = x / u
        y = y / u
        xmin = minimum(x); xmax = maximum(x)
        ymin = minimum(y); ymax = maximum(y)
    end

    myround(x; kwargs...) = oneunit(x) * round(x / oneunit(x); kwargs...)

    plot!([xmax, xmin, xmin, xmax, xmax],
        [ymax, ymax, ymin, ymin, ymax], color=:green)

    xticks = myround.([xmin, zero(xmin), xmax]; digits = 0)
#@show xmin, xmax, ymin, ymax
#@show xticks
    plot!(xtick = myround.([xmin, zero(xmin), xmax]; digits = 0))
    plot!(ytick = myround.([ymin, zero(xmin), ymax]; digits = 2))

    θ = LinRange(0, 2π, 501)
    rfov = sg.rfov
    scat!([0], [0], marker=:circle)
    plot!(rfov * cos.(θ), rfov * sin.(θ), color=:magenta) # rfov circle
    rfov = myround(sg.rfov, digits=1)
    title = "$(typeof(sg))"
    title = title[1:findfirst('{', title)-1]
    plot!(xlabel="x", ylabel="y"; title)

#=
    if sg isa SinoPar
    end
=#

    if sg isa SinoFan
        x0 = zero(sg.dso)
        y0 = sg.dso
        t = LinRange(0, 2π, 100)
        rot = sg.ar[1]
        rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
        p0 = rot * [x0; y0]
        pd = rot * [sg.xds'; sg.yds'] # detector points
        plot!(xlims = (-1,1) .* sg.dso ./ u)
        plot!(ylims = (-1,1) .* sg.dso ./ u)

        tmp = sg.ar .+ π/2 # trick: angle beta defined ccw from y axis
        scat!([p0[1]], [p0[2]], color=:blue, marker=:square) # source
        plot!(sg.dso * cos.(t), sg.dso * sin.(t), color=:cyan) # source circle
        scat!(sg.dso * cos.(tmp), sg.dso * sin.(tmp),
            color=:blue, marker=:circle, markersize=2) # source points
        scat!(vec(pd[1,:]), vec(pd[2,:]), marker=:circle,
            color=:orange, markersize=1) # detectors

        plot!([pd[1,1], p0[1], pd[1,end]], [pd[2,1], p0[2], pd[2,end]],
            color=:red, label="")
#       plot!(title="$(typeof(sg))")
    end

    if sg isa SinoMoj
        θ = LinRange(0, 2π, 100)
        rphi = sg.nb/2 * sg.d_moj.(θ)
        plot!(rphi .* cos.(θ), rphi .* sin.(θ), color=:blue, label="")
    #   rmax = maximum(sg.s)
    #   axis([-1 1 -1 1] * max([rmax ig.fov/2]) * 1.1)
    end


    return plot!(aspect_ratio=1)
end

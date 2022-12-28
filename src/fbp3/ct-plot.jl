#=
ct-plot.jl
Show CT geometries.
2022-10-01, Jeff Fessler
=#

export ct_geom_plot2!, ct_geom_plot3

#using Sinograms: CtGeom, RealU, source_zs
using .Plots: plot, plot!, default
using ImageGeoms: axes


function ct_geom_plot2!(st::CtGeom)
    plot!(; title = nameof(typeof(st)))
    if st isa CtFan
        sino_geom_plot_fan!(_ar(st), _rfov(st), _dso(st), _xds(st), _yds(st))
    end
    return plot!(aspect_ratio=1)
end


"""
    ct_geom_plot2!(st::CtGeom [, ig::ImageGeom])
Plot central 2D portion of the CBCT geometry.
"""
function ct_geom_plot2!(st::CtGeom, ig::ImageGeom)
    sino_geom_plot_ig!(_rfov(st); ig)
    ct_geom_plot2!(st)
end



function ct_geom_plot_ig!(ig::ImageGeom{3})
    default(label = "")
    xmin, xmax = extrema(axes(ig)[1])
    ymin, ymax = extrema(axes(ig)[2])
    zmin, zmax = extrema(axes(ig)[3])
    xticks = sort([zero(xmin), xmin, xmax])
    yticks = sort([zero(ymin), ymin, ymax])
    zticks = sort([zero(zmin), zmin, zmax])

    color = :green
    plot!(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin],
        ones(5) * zmin; color,
        xticks, yticks, zticks,
    )

    plot!(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin],
        ones(5) * zmax; color,
    )

    plot!(xmin*[1,1], ymin*[1,1], [zmin, zmax]; color)
    plot!(xmax*[1,1], ymin*[1,1], [zmin, zmax]; color)
    plot!(xmin*[1,1], ymax*[1,1], [zmin, zmax]; color)
    plot!(xmax*[1,1], ymax*[1,1], [zmin, zmax]; color)
end


"""
    ct_geom_plot3(st::CtFan [, ig::ImageGeom])
3D CBCT geometry plot
"""
function ct_geom_plot3(st::CtFan, ig::ImageGeom)

    default(label = "")
    plot()
    ct_geom_plot_ig!(ig)

    t1 = -_dso(st) * sin.(_ar(st))
    t2 = _dso(st) * cos.(_ar(st))
    t3 = _source_zs(st)
    plot!(t1, t2, t3, color = :blue, marker = :circle) # source points

    r100(x) = _round(maximum(abs, x); digits=2)
    r10(x) = _round(maximum(abs, x); digits=1)

    u = oneunit(_dso(st))
    plot!([-1, 1] * r100(t1), [0,0]*u, [0,0]*u, color=:red) # axes
    plot!([0,0]*u, [-1, 1]*r100(t2), [0,0]*u, color=:red)
    plot!([0,0]*u, [0,0]*u, [-1,1]*r10(t3), color=:red)

    detcolor = [:red, :cyan, :magenta]
    ia_list = unique(floor.(Int, 1 .+ (0:2)/3 .* st.na))

    for (ic, ia) in enumerate(ia_list)
        src = [t1[ia], t2[ia], t3[ia]]
        unit_vec = [-src[1], -src[2], 0u] / sqrt(src[1]^2 + src[2]^2)
        det_cen_loc = src + st.dsd * unit_vec
        plot!(map(x -> [x], det_cen_loc)..., marker=:star, color=:black)

        rot = _ar(st)[ia]
        rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
        pd = rot * [_xds(st)'; _yds(st)']

        color = detcolor[ic]
        for it in 1:st.nt
            plot!(pd[1,:], pd[2,:], src[3] .+ _t(st)[it] * ones(st.ns); color)
        end
    end
    plot!(aspect_ratio = :equal)
end

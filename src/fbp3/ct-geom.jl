#=
fbp3/ct-geom.jl
CT geometry for 3D image reconstruction
that describes the sampling characteristics of a CT scan
using the parallel, fan beam arc, or flat fan beam geometry.
Use sino_geom() instead for the 2D version.

todo: isolate the source trajectory? (e.g., helix, etc.)

2022-05-11, Jason Hu, Jeff Fessler translated from ct_geom.m in MIRT
=#

export CtGeom, CtParallel, CtFan
export CtPar, CtFanArc, CtFanFlat
export ct_geom_ge1
export dims

import Base: ones, zeros


"""
    CtGeom{Td,To}

Abstract type for representing 3D CT imaging geometries.
Currently supports only constant pitch source trajectories.

# Common struct elements
* `ns` size of each projection view
* `nt`
* `na` # of projection views
* `orbit` source orbit in degrees
* `orbit_start` starting angle in degrees
* `ds` detector pixel spacing
* `dt`
* `offset_s::Float32` # detector center offset (usually 0)
* `offset_t::Float32`

# Additional elements for `CtFan` types:
* `dsd` distance from source to detector
* `dso` distance from source to origin
* `dod` distance from origin to detector
* `pitch` helix pitch
* `source_z0` initial z position of x-ray source

# Basic methods
* `dims(cg)`
todo

#  Derived values (available by `getproperty`), i.e., `cg.?`

* `.s` s sample locations
* `.t` t sample locations
* `.gamma` [nb] gamma sample values [radians]
* `.gamma_s` gamma values given s values
* `.gamma_max` half of fan angle [radians], if offset_s=0
* `.gamma_max_abs` half of fan angle [radians] - recommended
* `.ws` (ns-1)/2 + st.offset_s
* `.wt` (nt-1)/2 + st.offset_t
* `.ad` [na] source angles in degrees
* `.ar` [na] source angles in radians
* `.dim` dimensions: [st.ns st.nt st.na]
* `.downsample(down)` reduce sampling by integer factor
* `.ones` ones(ns,nt,na, 'single')
* `.zeros` zeros(ns,nt,na, 'single')
* `.rad` [ns] radial distance of each ray from origin
* `.rfov` max radius within FOV
* `.footprint_size(ig)` max footprint width in 's' (not done)
* `.cone_angle` (half) cone angle on axis (s=0): +/- angle
* `.zfov` axial FOV
* `.source_zs` [na] z-locations of source for each view

# Methods
* `.down(down)` reduce sampling by integer factor
* `.shape(proj)` reshape `proj`
* `.plot2(ig)` plot the 2D geometry of the image
* `.plot3(ig)` plot the 3D geometry of the cone beam setup
* `.unitv(is,it,ia)` unit vector with single nonzero element
"""
abstract type CtGeom{Td,To} end

abstract type CtParallel{Td,To} <: CtGeom{Td,To} end #line 238 of ct_geom
abstract type CtFan{Td,To} <: CtGeom{Td,To} end


# types


struct CtPar{Td,To} <: CtParallel{Td,To}
    ns::Int
    nt::Int
    na::Int
    orbit::To
    orbit_start::To
    ds::Td
    dt::Td
    offset_s::Float32
    offset_t::Float32
end


struct CtFanArc{Td,To} <: CtFan{Td,To}
    ns::Int
    nt::Int
    na::Int
    orbit::To
    orbit_start::To
    ds::Td
    dt::Td
    offset_s::Float32
    offset_t::Float32
    dsd::Td
#   dso::Td
    dod::Td
#   dfs::Td
    pitch::Float32
#   user_source_zs::Vector{Td}
    source_z0::Td
end


struct CtFanFlat{Td,To} <: CtFan{Td,To}
    ns::Int
    nt::Int
    na::Int
    orbit::To
    orbit_start::To
    ds::Td
    dt::Td
    offset_s::Float32
    offset_t::Float32
    dsd::Td
#   dso::Td
    dod::Td
#   dfs::Td
    pitch::Float32
#   user_source_zs::Vector{Td}
    source_z0::Td
end


# constructors


function CtPar( ;
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
#   orbit::Union{Symbol,Real} = 360,
    orbit::RealU = 360,
    orbit_start::RealU = zero(orbit),
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    down::Int = 1,
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt)
    st = CtPar(ns, nt, na,
            To(orbit), To(orbit_start),
            Td(ds), Td(dt), Float32(offset_s), Float32(offset_t),
    )
    if down > 1
        st = downsample(st, down)
    end
    return st
end


function CtFanArc( ;
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
#   orbit::Union{Symbol,Real} = 360,
    orbit::RealU = 360,
    orbit_start::RealU = zero(orbit),
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    dsd::RealU = 4ns * oneunit(ds),
    dod::RealU = zero(ds),
    down::Int = 1,
    pitch::Real = 0,
    source_z0::RealU = zero(ds),
#   user_source_zs = [],
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt, dsd, dod, source_z0)
    st = CtFanArc(ns, nt, na,
            To(orbit), To(orbit_start),
            Td(ds), Td(dt), Float32(offset_s), Float32(offset_t), Td(dsd), Td(dod),
            Float32(pitch), Td(source_z0), # Td.(user_source_zs,
    )
    if down > 1
        st = downsample(st, down)
    end
    return st
end


function CtFanFlat( ;
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
#   orbit::Union{Symbol,Real} = 360,
    orbit::RealU = 360,
    orbit_start::RealU = zero(orbit),
    ds::RealU = 1,
    dt::RealU = ds,
    offset_s::Real = 0,
    offset_t::Real = 0,
    dsd::RealU = 4ns * oneunit(ds),
    dod::RealU = zero(ds),
    down::Int = 1,
    pitch::Real = 0,
    source_z0::RealU = zero(ds),
#   source_zs = (0:st.na-1) * source_dz .+ source_z0
)

    To = _promoter(orbit, orbit_start)
    Td = _promoter(ds, dt, dsd, dod, source_z0)
    st = CtFanFlat(ns, nt, na,
            To(orbit), To(orbit_start),
            Td(ds), Td(dt), Float32(offset_s), Float32(offset_t), Td(dsd), Td(dod),
            Float32(pitch), Td(source_z0), # Td.(user_source_zs,
    )
    if down > 1
        st = downsample(st, down)
    end
    return st
end



"""
    ct_geom_ge1( ; na::Int = 1)
CT geometry for GE Lightspeed system.
These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl.
"""
function ct_geom_ge1( ; na::Int = 1)
    return CtFanArc( ; ns = 888, nt = 64, na = 984, ds = 1.0239, dt = 1.0964,
        offset_s = 1.25, dsd = 949.075, dod = 408.075,
)
end


# methods that are common to all types

dims(st::CtGeom) = (st.ns, st.nt, st.na)
# angles() (in sino-geom.jl)
orbit_short(st::CtGeom) = 180 + 2 * rad2deg(st.gamma_max)
Base.ones(st::CtGeom) = ones(Float32, dims(st))::Array{Float32,3}
Base.zeros(st::CtGeom) = zeros(Float32, dims(st))::Array{Float32,3}

ct_geom_ws(st::CtGeom) = ((st.ns-1)//2 + st.offset_s)::Toffset
ct_geom_wt(st::CtGeom) = ((st.nt-1)//2 + st.offset_t)::Toffset
#ct_geom_s(st::CtGeom) = st.ds * ((0:st.ns-1) .- st.ws) # can't infer
function ct_geom_s(
    st::CtGeom{Td} ;
    T::DataType = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.ds, st.ws, st.ns)
end
function ct_geom_t(
    st::CtGeom{Td} ;
    T::DataType = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.dt, st.wt, st.nt)
end

# basic methods

# just to make the other parts of the code more convenient that use dfs
ct_geom_dfs(st::CtFanArc{Td}) where Td = zero(Td)
ct_geom_dfs(st::CtFanFlat{Td}) where Td = Inf * oneunit(Td)

ct_geom_rfov(st::CtPar) = maximum(abs, st.s)
ct_geom_rfov(st::CtFanArc) = st.dso * sin(ct_geom_gamma_max_abs(st))
ct_geom_rfov(st::CtFanFlat) = st.dso * sin(ct_geom_gamma_max_abs(st)) # todo: correct?

ct_geom_xds(st::CtFanArc) = st.dsd * sin.(st.gamma_s(st.s))
ct_geom_xds(st::CtFanFlat) = st.s
ct_geom_xds(st::CtParallel) = st.s

ct_geom_yds(st::CtFanArc) = @. st.dso - st.dsd * cos(st.gamma_s(st.s))
ct_geom_yds(st::CtFanFlat) = @. fill(-st.dod, st.ns)
ct_geom_yds(st::CtParallel) = 0 * st.s

ct_geom_cone_angle(st::CtParallel) = 0
ct_geom_cone_angle(st::CtFan) = atan((st.nt * st.dt)/2 / st.dsd)
#ct_geom_cone_angle(st::CtFan) = st.dso / st.dsd * st.nt * st.dt

ct_geom_zfov(st::CtParallel) = st.nt * st.dt
ct_geom_zfov(st::CtFan) = st.dso / st.dsd * st.nt * st.dt

ct_geom_gamma_s(st::CtFanArc, ss) = ss / st.dsd
ct_geom_gamma_s(st::CtFanFlat, ss) = atan.(ss / st.dsd)

ct_geom_gamma(st::CtGeom) = ct_geom_gamma_s(st, st.s)
ct_geom_gamma_max(st::CtGeom) = maximum(st.gamma)
ct_geom_gamma_max_abs(st::CtGeom) = maximum(abs, st.gamma)


# type inference help:
_rays_type(Td,To) = Iterators.ProductIterator{
    Tuple{LinRange{Td, Int}, LinRange{Td, Int}, LinRange{To, Int}, To}
}

"""
    i = rays(st::CtGeom)
Return parallel-beam coordinates of all rays for this CT geometry.
Return type of `i` is a `ProductIterator` that makes tuples of the form
`(u, v, ϕ, θ)`.
To make projections call
`p = [fun(c...) for c in i]` where `fun` is `radon(...)`.
"""
function rays(st::CtPar{Td,To})::_rays_type(Td,To) where {Td,To}
    u = st.s
    v = st.t
    ϕ = st.ar / oneunit(eltype(st.ar))
#   θ = [zero(To)]
    θ = zero(eltype(ϕ))
    i = Iterators.product(u, v, ϕ, θ)
    return i
#=
    u = [p[1] for p in i]
    v = [p[2] for p in i]
    ϕ = [p[3] for p in i]
    θ = [p[4] for p in i]
    Tϕ = eltype(oneunit(to_radians([oneunit(To)])[1]))
    return (u, v, ϕ, θ)::Tuple{Array{Td,4}, Array{Td,4}, Array{Tϕ,4}, Array{Tϕ,4}}
=#
end

#function cb_to_par(ss, tt, β)

function rays(st::CtFan{Td,To}) where {Td,To}
    st.pitch == 0 || throw("pitch not done")
    s = st.s
    t = st.t
    β = st.ar / oneunit(eltype(st.ar))
    i = Iterators.product(s, t, β)
    if st isa CtFanArc
        fun = stb -> cb_arc_to_par(stb..., st.dso, st.dod)
    else
        fun = stb -> cb_flat_to_par(stb..., st.dso, st.dod)
    end
    return Iterators.map(fun, i)
#=
    s = [p[1] for p in i]
    t = [p[2] for p in i]
    β = [p[3] for p in i]
    out = cb_arc_to_par.(s, t, st.dso, st.dod)
#   γ = ct_geom_gamma_s(st, s)
#   r = st.dso * sin.(γ)
#   ϕ = γ + β
#   θ =
    return (u, v, ϕ, θ)::Tuple{Array{Td,4}, Array{Td,4}, Array{Tϕ,4}, Array{Tϕ,4}}
=#
end


function ct_geom_source_dz_per_view(st::CtGeom{Td}) where Td
    if st.na == 1 || st.pitch == 0
        return zero(Td)
    end
#   if length(st.orbit) != 1 || st.orbit == 0
#       error("ERROR OCCURRED") # todo: cut?
#       return 0
#   end
    na_per_360 = st.na * (360 / st.orbit)
    out = st.pitch * st.zfov / na_per_360
end

function ct_geom_source_zs(st::CtGeom)
#   if !isempty(st.user_source_zs)
#       return st.user_source_zs
#   else
        source_dz = st.source_dz_per_view
        return (0:st.na-1) * source_dz .+ st.source_z0
#   end
end

#=
function ct_geom_footprint_size(st::CtPar, ig)
    di = sqrt(ig.dx^2 + ig.dy^2)
    smax = maximum(abs.(st.s))
    rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2
    dso = st.dso
    dsd = st.dsd
    if rfov > 0.99 * dso
        error("BAD DSO ERROR")
        return 1 / 0
    end
    return di / st.ds
end

function ct_geom_footprint_size(st::CtFanFlat, ig)
    di = sqrt(ig.dx^2 + ig.dy^2)
    smax = maximum(abs.(st.s))
    rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2
    dso = st.dso
    dsd = st.dsd
    if rfov > 0.99 * dso
        error("BAD DSO ERROR")
        return 1 / 0
    end
    return di / st.ds * sqrt(dsd^2 + smax^2) / (dso-rfov)
end

function ct_geom_footprint_size(st::CtFanArc, ig)
    di = sqrt(ig.dx^2 + ig.dy^2)
    smax = maximum(abs.(st.s))
    rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2
    dso = st.dso
    dsd = st.dsd
    if rfov > 0.99 * dso
        error("BAD DSO ERROR")
        return 1 / 0
    end
    return di / st.ds * dsd / (dso-rfov)
end
=#


function _downsample(st::CtGeom, down_s::Int, down_t::Int, down_a::Int)
    ns = 4 * max(st.ns ÷ 4down_s, 1)
    nt = 2 * max(st.nt ÷ 2down_t, 1)
    na = max(st.na ÷ down_t, 1)

#=
    user_spec = false
    if !isempty(st.user_source_zs)
        st.user_source_zs = st.user_source_zs[1:down_a:st.na]
        user_spec = true
    end

    if length(st.orbit_start) > 1
        st.orbit_start = st.orbit_start[1:down_a:si.na]
        user_spec = true
    end
=#

#=
    if user_spec
        st.na = length([1:down_a:st.na])
    elseif all(diff(si.source_zs) == 0)
        st.na = max(1, round(si.na / down_a))
    else
        nturn = si.orbit / 360
        na1 = si.na / nturn

        tol = 0.0001
        if abs(round(na1)-na1) < tol
            na2 = round(na1 / down_a)
            tmp = nturn * na2
            if abs(round(tmp)-tmp) < tol
                st.na = round(tmp)
            else
                st.na = floor(tmp)
                st.orbit = 360 / na2 * st.na
            end
        else
            st.na = round(si.na / down_a)
        end
    end
=#

    out = (ns, nt, na, st.orbit, st.orbit_start,
         st.ds * down_s, st.dt * down_t, st.offset_s, st.offset_t)
    if st isa CtFan
         out = (out..., st.dsd, st.dod, st.pitch, st.source_z0)
    end
    return out
end


"""
    downsample(ct, down::Int)
    downsample(ct, down::NTuple{3,Int})
Down-sample CT geometry (for testing with small problems).
"""
downsample(ct::CtGeom, down::Int) = downsample(ct, (down,down,down))

function downsample(ct::C, down::NTuple{3,Real}) where {C <: CtGeom}
    return all(==(1), down) ? ct : C(_downsample(ct, down...)...)::C
end

function ct_geom_unitv(st::CtGeom ;
    is::Int = floor(Int, st.ns/2) + 1,
    it::Int = floor(Int, st.nt/2) + 1,
    ia::Int = 1,
)
    out = zeros(st)
    out[is, it, ia] = 1
    return out
end


#=
# plot 2D geometry of the setup

using PyPlot

function ct_geom_plot2(st::CtGeom, ig::ImageGeom)
    t = LinRange(0, 2*pi, 1001)

    if isa(st, CtParallel)
        error("Not implemented yet")
    end
    beta0 = deg2rad(st.orbit_start)
    rot = [cos(beta0) sin(beta0); -sin(beta0) cos(beta0)]
    p0 = st.dso * [-sin(beta0); cos(beta0)]
    pd = rot' * [st.xds'; st.yds']
    rfov = st.dso * sin(maximum(abs.(st.gamma)))

    Plots.scatter([p0[1]],[p0[2]], linestyle = :solid)
    Plots.plot!(st.dso * cos.(t), st.dso * sin.(t), linecolor = :red)
    Plots.plot!(pd[1,:], pd[2,:], linestyle = :solid, linecolor = :yellow)

    xmin = minimum(ig.x)
    xmax = maximum(ig.x)
    ymin = minimum(ig.y)
    ymax = maximum(ig.y)
    Plots.plot!([xmax, xmin, xmin, xmax, xmax], [ymax, ymax, ymin, ymin, ymax], linecolor = :green)

    Plots.plot!([pd[1,1], p0[1], last(pd[1,:])], [pd[2,1], p0[2], last(pd[2,:])], linestyle = :dash)
    Plots.plot!(rfov*cos.(t), rfov*sin.(t), linecolor = :magenta, linestyle = :dot)
end

#helper functions for ct_geom_plot3
function r100(x)
    return 100 * ceil(maximum(abs.(x))/100)
end

function r10(x)
    return 10 * ceil(maximum(abs.(x))/10)
end

#3D helical geometry plot, saves plot to myfig.png
function ct_geom_plot3(st::CtGeom, ig::ImageGeom)
    if isinf(st.dso) || isinf(st.dsd)
        error("Parallel beam Not done")
    end

    t1 = -st.dso * sin.(st.ar)
    t2 = st.dso * cos.(st.ar)
    t3 = st.source_zs

    fig = figure()
    ax = Axes3D(fig)
    ax.plot3D(t1, t2, t3, color = "yellow", linestyle = "dotted")

    ax.plot3D([-1, 1] * r100(t1), [0,0], [0,0])
    ax.plot3D([0,0], [-1, 1]*r100(t2), [0,0])
    ax.plot3D([0,0], [0,0], [-1,1]*r10(t3))

    #skip the part about ztick line 740
    xmin = minimum(ig.x)
    xmax = maximum(ig.x)
    ymin = minimum(ig.y)
    ymax = maximum(ig.y)
    zmin = minimum(ig.z)
    zmax = maximum(ig.z)
    ax.plot3D([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], [zmin, zmin, zmin, zmin, zmin], color = "green", linestyle = "dashed")
    ax.plot3D([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], [zmax, zmax, zmax, zmax, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmin*[1,1], ymin*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmax*[1,1], ymin*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmin*[1,1], ymax*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmax*[1,1], ymax*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmin*ones(ig.nz), ymin*ones(ig.nz), ig.z, color = "green", linestyle = "dotted")

    detcolor = ["red", "cyan", "magenta"]
    if st.na == 1
        ia_list = 1
    elseif st.na == 2
        ia_list = [1, st.na]
    else
        ia_list = [1, 1+floor(Int, st.na/2), st.na]
    end

    for ia = ia_list
        src = [t1[ia], t2[ia], t3[ia]]
        unit_vec = [-src[1], -src[2], 0]/sqrt(src[1]^2 + src[2]^2)
        det_cen_loc = src + st.dsd * unit_vec
        #ax.plot3D(det_cen_loc[1], det_cen_loc[2], det_cen_loc[3], )

        rot = st.ar[ia]
        rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
        pd = rot * [st.xds'; st.yds']

        for it = 1:st.nt
            indices = findall(x -> x == ia, ia_list)
            indices = indices[1]
            ax.plot3D(pd[1,:], pd[2,:], src[3] .+ st.t[it] * ones(st.ns), color = detcolor[indices])
        end
    end
    ax.view_init(22, -200)
    PyPlot.savefig("myfig.png")
end
=#


ct_geom_fun0 = Dict([
#   (:help, st -> ct_geom_help()),

#   (:dim, st -> (st.ns, st.nt, st.na)),
    (:ws, st -> ct_geom_ws(st)),
    (:wt, st -> ct_geom_wt(st)),
#   (:ones, st -> ones(Float64, st.dim)),
#   (:zeros, st -> zeros(Float64, st.dim)),
    (:unitv, st -> ((; kwarg...) -> ct_geom_unitv(st ; kwarg...))),

    (:s, st -> ct_geom_s(st)),
    (:t, st -> ct_geom_t(st)),
    (:dso, st -> st.dsd - st.dod),
    (:dfs, st -> ct_geom_dfs(st)),

    (:gamma, st -> ct_geom_gamma(st)),
    (:gamma_s, st -> (ss -> ct_geom_gamma_s(st, ss))),
    (:gamma_max, st -> ct_geom_gamma_max(st)),
    (:gamma_max_abs, st -> ct_geom_gamma_max_abs(st)),
    (:xds, st -> ct_geom_xds(st)),
    (:yds, st -> ct_geom_yds(st)),

    (:cone_angle, st -> ct_geom_cone_angle(st)),
    (:zfov, st -> ct_geom_zfov(st)),
    (:source_dz_per_view, st -> ct_geom_source_dz_per_view(st)),
    (:source_zs, st -> ct_geom_source_zs(st)),
    (:orbit_short, st -> orbit_short(st)),
#   (:downsample, st -> (down::Int -> ct_geom_downsample(st, down))),

    (:ad, st -> angles(st)),
    (:ar, st -> to_radians(st.ad)),
#   (:rad, st -> ct_geom_rad(st)),
    (:rfov, st -> ct_geom_rfov(st)),
    (:shape, st -> (proj -> reshaper(proj, dims(st)))),
#   (:footprint_size, st -> ct_geom_footprint_size(st)),

#   (:plot_grid, st -> ((plot::Function) -> ct_geom_plot_grid(st, plot))),
#   (:plot2, st -> ((ig::ImageGeom) -> ct_geom_plot2(st, ig))),
#   (:plot3, st -> ((ig::ImageGeom) -> ct_geom_plot3(st, ig)))
])

Base.getproperty(st::CtGeom, s::Symbol) =
    haskey(ct_geom_fun0, s) ? ct_geom_fun0[s](st) :
    getfield(st, s)

Base.propertynames(st::CtGeom) =
    (fieldnames(typeof(st))..., keys(ct_geom_fun0)...)


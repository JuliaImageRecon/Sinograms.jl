#=
fbp3/cbct-back.jl
Translated from cbct_back.m in MIRT
Copyright 2022-5-18 Jason Hu and Jeff Fessler, University of Michigan
=#

using FFTW
using ImageGeoms: ImageGeom
export cbct_back


"""
    cbct_back(proj, cg, ig)

Cone-beam backprojector for feldkamp.jl

# in
* `proj` (ns,nt,na)     cone-beam projection views
* `cg::CtGeom`
* `ig::ImageGeom`

# out
* `img` (nx,ny,nz) back-projection result
"""
function cbct_back(
    proj::AbstractArray{Ts,3},
    cg::CtFan{Td, To},
    ig::ImageGeom{3} ;
    ia_skip::Int = 1,
) where {Ts <:Number, Td, To}

    # type inference help:
    Toffset = Float32 # eltype(cg.offset_s)
    T = eltype(oneunit(Ts) *
        (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

    return cbct_back_fan(proj,
        cg.ar,
        cg.dsd, cg.dso,
#       cg.offset_source::RealU,
        cg.ds, cg.dt,
        cg.offset_s, cg.offset_t,
        cg isa CtFanArc, # is_arc
#       source_zs = zeros(na),
        axes(ig)...,
        ig.mask ;
        ia_skip,
    )::Array{T,3}
end

#=
    ns = cg.ns
    nt = cg.nt
    na = cg.na
    ds = cg.ds
    dt = cg.dt
    offset_s = cg.offset_s
    offset_t = cg.offset_t
    offset_source = 0
    dsd = cg.dsd
    dso = cg.dso
    dfs = cg.dfs
    orbit = cg.orbit
    orbit_start = cg.orbit_start
    source_zs = zeros(na) # in matlab this array is all zeros
    if isa(ig, im_geom2)
        mask = ig.mask
    elseif isa(ig, ImageGeom)
        mask = ig.mask_or
    else
        error("Mask Not implemented yet")
    end
    nz = ig.nz
    dx = ig.dx
    dy = ig.dy
    dz = ig.dz
    offset_xyz = [ig.offset_x, ig.offset_y, ig.offset_z]
    ia_skip = 1
    scale_dang = true

    (nx, ny) = size(mask)
    betas = @. deg2rad(orbit * (0:na-1) / na + orbit_start)
    wx = (nx-1)/2 + offset_xyz[1]
    wy = (ny-1)/2 + offset_xyz[2]
    wz = (nz-1)/2 + offset_xyz[3]

    xc, yc = ndgrid(((0:nx-1) .- wx) * dx, ((0:ny-1) .- wy) * dy)
    zc = ((0:nz-1) .- wz) * dz

    xc = xc[mask]
    yc = yc[mask]

    ws = (ns+1)/2 + offset_s
    wt = (nt+1)/2 + offset_t

    img = zeros(nx, ny, nz)
    sdim = (ns+1, nt)
    proj1 = zeros(ns+1, nt)

    for iz in 1:nz
        ia_min = 1
        ia_max = na

        img2 = zeros(ns*nt-1)
        for ia in ia_min:ia_skip:ia_max
            beta = betas[ia]

            x_beta = xc * cos.(beta) + yc * sin.(beta)
            y_betas = - (-xc * sin.(beta) + yc * cos.(beta)) .+ dso

            if isinf(dsd) || isinf(dso)
                mag = 0*y_betas .+ 1
            else
                mag = dsd ./ y_betas
            end

            if isinf(dfs) || isinf(dsd) || isinf(dso)
                sprime = mag .* x_beta
            elseif dfs == 0
                r_loop = x_beta .- offset_source
                sprime = dsd * atan.(r_loop, y_betas)
            end

            tprime = mag * (zc[iz] - source_zs[ia])

            bs = sprime / ds .+ ws
            bt = tprime / dt .+ wt

            bs[bs .< 1 .|| bs .> ns] .= ns+1
            bt = max.(bt, 1)
            bt = min.(bt, nt)

            is = floor.(bs)
            it = floor.(bt)

            is[is .== ns+1] .= ns
            it[it .== nt] .= nt - 1

            wr = bs - is
            wl = -wr .+ 1
            wu = bt - it
            wd = -wu .+ 1

            proj1[1:ns, :] = proj[:,:,ia]
            is = convert(Vector{Int64}, is)
            it = convert(Vector{Int64}, it)
            p1 = wl .* proj1[Base._sub2ind(sdim, is, it)] + wr .* proj1[Base._sub2ind(sdim, is .+ 1, it)]
            p2 = wl .* proj1[Base._sub2ind(sdim, is, it .+ 1)] + wr .* proj1[Base._sub2ind(sdim, is .+ 1, it .+ 1)]

            p0 = wd .* p1 + wu .* p2

            if isinf(dfs) || isinf(dsd) || isinf(dso)
                p0 = p0 .* mag .^2
            elseif dfs == 0
                p0 = p0 .* (dsd.^2) ./ (r_loop.^2 + y_betas.^2)
            end

            img2 = img2 + p0
        end

        img[:,:,iz] = embed(img2, mask)
    end

    if scale_dang
        img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img
    end

    return img
end
=#


function cbct_back_fan(
    proj::AbstractArray{<:Ts,3},
    betas::AbstractVector{To},
    dsd::RealU,
    dso::RealU,
#   offset_source::RealU,
    ds::Tds,
    dt::Tds,
    offset_s::Toffset,
    offset_t::Toffset,
    is_arc::Bool,
#   source_zs = zeros(na),
    xc::AbstractVector{<:Tc},
    yc::AbstractVector{<:Tc},
    zc::AbstractVector{<:Tc},
    mask::AbstractArray{Bool,3} ;
    ia_skip::Int = 1,
    T::DataType = eltype(oneunit(Ts)
        * (oneunit(To) * oneunit(Tc) / oneunit(Tds) + oneunit(Toffset))),
)::Array{T,3} where {
    Ts <: Number,
    To <: RealU,
    Tds <: RealU,
    Toffset <: Real,
    Tc <: RealU,
}

    image = zeros(T, size(mask)) # need zero(T) outside mask

    cbct_back_fan!(image, proj, betas,
        dsd, dso,
#       offset_source,
        ds, dt, offset_s, offset_t, is_arc,
#       source_zs = zeros(na),
        xc, yc, zc, mask ; ia_skip,
    )

    return image
end


"""
    cbct_back_fan!(...)
This should work even for non-uniformly spaced source angles
(for a circular source trajectory).
"""
function cbct_back_fan!(
    image::Array{T,3},
    proj::AbstractArray{<:Number,3},
    betas::AbstractVector{To},
    dsd::RealU,
    dso::RealU,
#   offset_source::RealU,
    ds::RealU,
    dt::RealU,
    offset_s::Toffset,
    offset_t::Toffset,
    is_arc::Bool,
#   source_zs = zeros(na),
    xc::AbstractVector{<:RealU},
    yc::AbstractVector{<:RealU},
    zc::AbstractVector{<:RealU},
    mask::AbstractArray{Bool,3} ;
    ia_skip::Int = 1,
) where {T <: Number, To <: RealU, Toffset <: Real}

    length.((xc,yc,zc)) == size(image) == size(mask) || throw("size mismatch")

    (ns, nt, na) = size(proj)
    ws = Toffset((ns+1)/2 + offset_s)
    wt = Toffset((nt+1)/2 + offset_t)
    if ia_skip > 1
        ia = 1:ia_skip:na
        proj = @view proj[:,:,ia]
        betas = @view betas[ia]
    end
    sinβ = sin.(betas)
    cosβ = cos.(betas)

    # scale projections by dβ for Riemann-like integration
    dβ = diff(betas)
    dβ = [dβ[1]; dβ] # todo
    proj = proj .* reshape(dβ, 1, 1, :) / 2

    xc_ds = xc / ds
    yc_ds = yc / ds
    zc_ds = zc / ds

    image[.! mask] .= zero(T)
    Threads.@threads for c in findall(mask)
        image[c] = cbct_back_fan_voxel(
            proj, sinβ, cosβ, dt / ds, ws, wt,
#           offset_source / ds,
            dsd / ds, dso / ds, is_arc,
            xc_ds[c[1]], yc_ds[c[2]], zc_ds[c[3]] ;
            T
        )
    end

    return image
end


# back-project all views to one voxel
# at location (xc,yc,zc)
function cbct_back_fan_voxel(
    proj::AbstractArray{Tp,3},
    sinβ::AbstractVector{To},
    cosβ::AbstractVector{To},
    dt_ds::Real,
    ws::Tw,
    wt::Tw,
#   offset_source_ds::Real,
    dsd_ds::RealU,
    dso_ds::RealU,
    is_arc::Bool,
#   source_zs = zeros(na) #for some reason in matlab this array is all zeros
    x_ds::Tx,
    y_ds::Tx,
    z_ds::Tx ;
    T::DataType = eltype(oneunit(Tp) * one(To) * one(Tw) * one(Tx)),
) where {Tp <: Number, To <: Real, Tw <: Real, Tx <: Real}

    voxel = zero(Tp)
    (ns, nt, na) = size(proj)

    for ia in 1:na
        @inbounds sβ = sinβ[ia]
        @inbounds cβ = cosβ[ia]

        x_beta = x_ds * cβ + y_ds * sβ
        y_beta = dso_ds - (-x_ds * sβ + y_ds * cβ) # dso - y_beta

        if is_arc
            r_loop = x_beta # - offset_source_ds
            sprime_ds = dsd_ds * atan(r_loop, y_beta)
            mag = dsd_ds / y_beta
            w2 = abs2(dsd_ds) / (abs2(r_loop) + abs2(y_beta))
        else
            mag = dsd_ds / y_beta
            sprime_ds = mag * x_beta
            w2 = abs2(mag) # image weighting
        end

#       tprime_ds = mag * (z_ds - source_zs[ia])
        tprime_ds = mag * (z_ds - 0)

        bs = sprime_ds + ws # unitless bin index
        bt = tprime_ds / dt_ds + wt

        is = floor(Int, bs)
        it = floor(Int, bt)

        if (1 ≤ is < ns) && (1 ≤ it < nt)

            # bilinear interpolation:
            wr = bs - is # right weight
            wl = 1 - wr # left weight
            wu = bt - it
            wd = 1 - wu

            @inbounds p1 = wl * proj[is, it, ia] + wr * proj[is + 1, it, ia]
            @inbounds p2 = wl * proj[is, it + 1, ia] + wr * proj[is + 1, it + 1, ia]

            voxel += w2 * (wd * p1 + wu * p2)
        end
    end

    return voxel
end

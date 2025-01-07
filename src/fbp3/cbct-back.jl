#=
fbp3/cbct-back.jl
Translated from cbct_back.m in MIRT
Copyright 2022-5-18 Jason Hu and Jeff Fessler, University of Michigan
=#

using FFTW
using ImageGeoms: ImageGeom
export cbct_back


"""
    cbct_back(proj, rg, ig)

Cone-beam backprojector for feldkamp.jl

# in
* `proj (ns,nt,na)` cone-beam projection views
* `rg::CtGeom`
* `ig::ImageGeom`

# out
* `img (nx,ny,nz)` back-projection result
"""
function cbct_back(
    proj::AbstractArray{<:Number,3},
    rg::CtFan,
    ig::ImageGeom{3},
    ;
    ia_skip::Int = 1,
)

    return cbct_back_fan(proj,
        _ar(rg), # "betas"
        rg.dsd, _dso(rg),
#       rg.offset_source::RealU,
        rg.ds, rg.dt,
        rg.offset_s, rg.offset_t,
        rg isa CtFanArc, # is_arc
#       source_zs = zeros(na),
        axes(ig)...,
        ig.mask ;
        ia_skip,
    )
end


function cbct_back_fan(
    proj::AbstractArray{Ts,3},
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
    xc::AbstractVector{Tc},
    yc::AbstractVector{Tc},
    zc::AbstractVector{Tc},
    mask::AbstractArray{Bool,3},
    ;
    ia_skip::Int = 1,
    T::Type{<:Number} = typeof(oneunit(Ts) *
        (oneunit(To) * oneunit(Tc) / oneunit(Tds) + oneunit(Toffset))),
) where {
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
    betas::AbstractVector{<:RealU},
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
    mask::AbstractArray{Bool,3},
    ;
    ia_skip::Int = 1,
) where {T <: Number, Toffset <: Real}

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

    xc_ds = xc / ds
    yc_ds = yc / ds
    zc_ds = zc / ds

    image[.! mask] .= zero(T)
    Threads.@threads for c in findall(mask)
        image[c] = cbct_back_fan_voxel(
            proj, sinβ, cosβ, dt / ds, ws, wt,
#           offset_source / ds,
            dsd / ds, dso / ds, is_arc,
            xc_ds[c[1]], yc_ds[c[2]], zc_ds[c[3]],
        )
    end # COV_EXCL_LINE

    return image
end


# back-project all views to one voxel
# at location (xc,yc,zc)
# Output voxel has same units as `proj`.
@inline function cbct_back_fan_voxel(
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
    z_ds::Tx,
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

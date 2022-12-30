# fbp2/back-fan.jl

export fbp_back

using LazyGrids: ndgrid
using ImageGeoms: ImageGeom, embed
#using Sinograms: SinoFan


# fan-beam case
function fbp_back(
    rg::SinoFan{Td, To},
    ig::ImageGeom,
    sino::AbstractMatrix{Ts} ;
    ia_skip::Int = 1,
) where {Td, To, Ts <: Number}

    dims(rg) == size(sino) || throw("sino size")

    is_arc = iszero(_dfs(rg)) ? true : isinf(_dfs(rg)) ? false : throw("bad dfs")

    # type inference help:
    Toffset = Float32 # typeof(rg.offset)
    T = typeof(oneunit(Ts) * (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

    return fbp_back_fan(
        sino, _ar(rg),
        rg.dsd, _dso(rg),
        rg.source_offset, is_arc,
        rg.d, rg.offset,
#       _rfov(rg),
#       ndgrid(axes(ig)...)...,
        axes(ig)...,
        ig.mask ; ia_skip,
    )::Matrix{T}
end


#=
# old "matlab-like" way with lots of broadcast
function fbp_back_fan_old(
    sino::AbstractMatrix{Ts},
    betas::AbstractVector{To},
    dsd::RealU,
    dso::RealU,
    source_offset::RealU,
    is_arc::Bool,
    ds::Tds,
    offset::Toffset,
#   rfov::RealU,
    xc::AbstractArray{Tc},
    yc::AbstractArray{Tc},
    mask::AbstractMatrix{Bool} ;
    ia_skip::Int = 1,
    warned::Bool = false,
) where {Tds <: RealU, Tc <: RealU, Toffset <: Real, Ts <: Number, To <: RealU}

    Td = promote_type(Tds, Tc)
    T = typeof(oneunit(Ts) * (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

    nb, na = size(sino)

    # trick: extra zero column saves linear interpolation indexing within loop!
    sino = cat(dims=1, sino, zeros(Ts, 1, size(sino,2)))

#=
    if do_r_mask
        rr = @. sqrt(abs2(xc) + abs2(yc)) # (nx,ny)
        mask = mask .& (rr .< rfov)
    end
=#

    xc = xc[vec(mask)] # [np] pixels within mask
    yc = yc[vec(mask)]

    wb = Toffset((nb+1)/2 + offset)

    img = zero(T)

    for ia in 1:ia_skip:na

        beta = betas[ia]
        d_loop = @. (dso + xc * sin(beta) - yc * cos(beta)) # dso - y_beta
        r_loop = @. (xc * cos(beta) + yc * sin(beta) - source_offset) # x_beta-roff

        if is_arc
            sprime_ds = (dsd/ds) * atan.(r_loop, d_loop) # s' / ds
            w2 = @. (abs2(dsd) / (abs2(d_loop) + abs2(r_loop))) # [np] image weighting
        else # flat
            mag = dsd ./ d_loop
            sprime_ds = mag .* r_loop / ds
            w2 = abs2.(mag) # [np] image-domain weighting
        end

        bb = sprime_ds .+ wb # [np] bin "index"

#=
        # nearest neighbor interpolation:
        ib = round.(Int, bb)
        # trick: make out-of-sinogram indices point to those extra zeros
        @. ib[ib < 1 | ib > nb] = nb+1
        img += sino[ib, ia] ./ w2
=#

        # linear interpolation:
        il = floor.(Int, bb) # left bin
        ir = 1 .+ il # right bin

        # handle truncated sinograms using the extra column of zeros
        ig = (il .≥ 1) .& (ir .≤ nb)
        if !all(ig)
            if !warned
                @warn("image exceeds system FOV; modify ig.mask?")
                warned = true
            end
            il[.!ig] .= nb+1
            ir[.!ig] .= nb+1
        end

        wr = bb - il # left weight
        wl = 1 .- wr # right weight

        img = @. (img + (wl * sino[il, ia] + wr * sino[ir, ia]) * w2)
    end

    img .*= (π * ia_skip / na)
    return embed(img, mask)::Matrix{T}
end
=#


# new version with threads and simplified code

"""
    fbp_back_fan(sino, betas,
        dsd, dso, dfs, source_offset, is_arc,
        ds, offset, xc, yc, mask ; ia_skip, T)
Pixel-driven back-projection
for a grid of `(xc,yc)` pixel center locations
for sinogram `sino` from a
fan-beam geometry.
It assumes the angles are equally spaced over `[0,π)`.

# in
- `sino::Matrix{<:Number}` `(nb, na)` usually ramp-filtered
- `betas::Vector{<:Real}` `(na)` in radians
- `dsd,dso,source_offset::RealU` geometry
- `is_arc::Bool` arc or flat?
- `ds::RealU` ray spacing
- `offset::Real` detector offset (usually 0)
- `xc::Vector{<:RealU}` `(nx)` pixel centers
- `yc::Vector{<:RealU}` `(ny)` pixel centers
- `mask::Matrix{Bool}` `(nx, ny)` which pixels to reconstruct

# option
- `ia_skip::Int` default 1
- `T::Type{<:Number}` usually same as `eltype(sino)`

# out
- `image::Matrix` `(nx, ny)`
"""
function fbp_back_fan(
    sino::AbstractMatrix{Ts},
    betas::AbstractVector{To},
    dsd::RealU,
    dso::RealU,
    source_offset::RealU,
    is_arc::Bool,
    ds::Tds,
    offset::Toffset,
    xc::AbstractArray{Tc},
    yc::AbstractArray{Tc},
    mask::AbstractMatrix{Bool} ;
    ia_skip::Int = 1,
    T::Type{<:Number} = typeof(oneunit(Ts) *
        (oneunit(To) * oneunit(Tc) / oneunit(Tds) + oneunit(Toffset))),
) where {Ts <: Number, To <: RealU, Tds <: RealU, Toffset <: Real, Tc <: RealU}

    image = zeros(T, size(mask)) # need zero(T) outside mask
    sinβ = sin.(betas[1:ia_skip:end])
    cosβ = cos.(betas[1:ia_skip:end])
    fbp_back_fan!(image, sino, sinβ, cosβ,
        dsd, dso, source_offset, is_arc,
        ds, offset, xc, yc, mask ; ia_skip)
    return image
end


"""
    fbp_back_fan!(image, sino, sinβ, cosβ,
        dsd, dso, source_offset, is_arc,
        ds, offset, xc, yc, mask ; ia_skip)
Mutating version of
pixel-driven back-projection
for a grid of `(xc,yc)` pixel center locations
for sinogram `sino` from a fan-beam geometry.
It uses `Threads`.
It assumes the angles are equally spaced over `[0,π)`.

# in
- `sino::Matrix{<:Number}` `(nb, na)` usually ramp-filtered
- `sinβ::Vector{<:Real}` `(na)`
- `cosβ::Vector{<:Real}` `(na)`
- `dsd,dso,source_offset::RealU` geometry
- `is_arc::Bool`
- `ds::RealU` ray spacing
- `offset::Real` detector offset (usually 0)
- `xc::Vector{<:RealU}` `(nx)` pixel centers
- `yc::Vector{<:RealU}` `(ny)` pixel centers
- `mask::Matrix{Bool}` `(nx, ny)` which pixels to reconstruct

# option
- `ia_skip::Int` default 1

# out
- `image::Matrix` `(nx, ny)` matrix to be mutated
"""
function fbp_back_fan!(
    image::AbstractMatrix{T},
    sino::AbstractMatrix{<:Number},
    sinβ::AbstractVector{<:Real}, # sin.(angles) (na_subset)
    cosβ::AbstractVector{<:Real}, # cos.(angles) (na_subset)
    dsd::RealU,
    dso::RealU,
    source_offset::RealU,
    is_arc::Bool,
    ds::RealU,
    offset::Toffset,
    xc::AbstractArray{<:RealU},
    yc::AbstractArray{<:RealU},
    mask::AbstractMatrix{Bool} ;
    ia_skip::Int = 1,
) where {T <: Number, Toffset <: Real}

    length.((xc,yc)) == size(image) == size(mask) || throw("size mismatch")

    nb = size(sino, 1)

    wb = Toffset((nb+1)/2 + offset)
    if ia_skip > 1
        sino = @view sino[:,1:ia_skip:end]
    end
    xc_ds = xc / ds
    yc_ds = yc / ds

    image[.! mask] .= zero(T)

    Threads.@threads for c in findall(mask)
        image[c] = fbp_back_fan_xy(
            sino, sinβ, cosβ,
            dsd/ds, dso/ds, source_offset/ds, is_arc,
            wb, xc_ds[c[1]], yc_ds[c[2]]; T,
        )
    end

    return image
end


"""
    fbp_back_fan_xy(sino, sinβ, cosβ,
        dsd_ds, dso_ds, source_offset_ds, is_arc,
        wb, x, y ; T)
Pixel-driven back-projection for a single (x,y) location
for sinogram `sino` from a
fan-beam geometry.
It assumes the angles are equally spaced over `[0,π)`.

# in
- `sino::Matrix{<:Number}` `(nb, na)` usually ramp-filtered
- `sinβ::Vector{<:Real}` `(na)`
- `cosβ::Vector{<:Real}` `(na)`
- `dsd_ds,dso_ds,source_offset_ds::Real` geometry, normalized
- `wb::Real = (nb+1)/2 + offset` where usually `offset=0`
- `x,y::Real` pixel center location, normalized by ray spacing

# option
- `T::Type{<:Number}` typically same as `eltype(sino)`

# out
- Returns a scalar of type `T`.
"""
function fbp_back_fan_xy(
    sino::AbstractMatrix{Ts}, # (nb, na_subset)
    sinβ::AbstractVector{To}, # sin.(angles) (na_subset)
    cosβ::AbstractVector{To}, # cos.(angles) (na_subset)
    dsd_ds::Real,
    dso_ds::Real,
    source_offset_ds::Real,
    is_arc::Bool,
    wb::Tb, # (nb+1)/2 + offset
    x_ds::Tx, # xc / ds
    y_ds::Tx ;
    T::Type{<:Number} = typeof(oneunit(Ts) * one(To) * one(Tb) * one(Tx)),
) where {Ts <: Number, To <: Real, Tb <: Real, Tx <: Real}

    nb = size(sino,1)
    na_subset = size(sino,2)

    pixel = zero(T)

    for ia in 1:na_subset
        @inbounds sβ = sinβ[ia]
        @inbounds cβ = cosβ[ia]
        d_loop = dso_ds + x_ds * sβ - y_ds * cβ # dso - y_beta
        r_loop = (x_ds * cβ + y_ds * sβ - source_offset_ds) # x_beta-roff

        if is_arc
            sprime_ds = dsd_ds * atan(r_loop, d_loop) # s' / ds
            w2 = abs2(dsd_ds) / (abs2(d_loop) + abs2(r_loop)) # image weighting
        else # flat
            mag = dsd_ds / d_loop
            sprime_ds = mag * r_loop # / ds
            w2 = abs2(mag) # image weighting
        end

        bb = sprime_ds + wb # unitless bin "index"

#=
        # nearest neighbor interpolation:
        ib = round(Int, bb)
        if 1 ≤ ib ≤ nb
            @inbounds pixel += sino[ib, ia] / w2
        end
=#

        # linear interpolation:
        il = floor(Int, bb) # left bin
        ir = 1 + il # right bin

        if (1 ≤ il) && (ir ≤ nb)
            wr = bb - il # right weight
            wl = 1 - wr # left weight
            @inbounds pixel += (wl * sino[il, ia] + wr * sino[ir, ia]) * w2
        end
    end

    return pixel # * (π / na_subset) todo
end

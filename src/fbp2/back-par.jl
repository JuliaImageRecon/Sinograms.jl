# fbp2/back-par.jl

export fbp_back

using LazyGrids: ndgrid
using ImageGeoms: ImageGeom, embed
#using Sinograms: SinoPar



"""
    img = fbp_back(rg, ig, sino ; ia_skip)

2D pixel-driven backprojection for FBP.

# in
- `rg::SinoGeom`
- `ig::ImageGeom`
- `sino::AbstractArray{<:Number}` sinogram(s) (line integrals), usually ramp filtered

# options
- `ia_skip::Int` downsample in angle to save time for quick tests (default: 1)

# out
- `img::Array{<:Number}` reconstructed image(s)
"""
fbp_back

# parallel-beam case
function fbp_back(
    rg::SinoPar{Td, To},
    ig::ImageGeom,
    sino::AbstractMatrix{Ts} ;
    ia_skip::Int = 1,
#   do_r_mask::Bool = false
) where {Td, To, Ts <: Number}

    dims(rg) == size(sino) || throw("sino size")

    # type inference help:
    Toffset = Float32 # typeof(rg.offset)
    T = typeof(oneunit(Ts) * (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

    return fbp_back_par(
        sino, _ar(rg),
        rg.d, rg.offset,
#       _rfov(rg),
#       ndgrid(axes(ig)...)...,
        axes(ig)...,
        ig.mask ; ia_skip,
    )::Matrix{T}
end


#=
# old "matlab-like" way with lots of broadcast
function fbp_back_par_old(
    sino::AbstractMatrix{Ts},
    angles::AbstractVector{To},
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

        angle = angles[ia]
        (sinϕ, cosϕ) = sincos(angle)
        bb = @. (xc * cosϕ + yc * sinϕ) # [np]
        bb = @. (bb / ds + wb) # unitless bin index

#=
        # nearest neighbor interpolation:
        ib = round.(Int, bb)
        # trick: make out-of-sinogram indices point to those extra zeros
        @. ib[ib < 1 | ib > nb] = nb+1
        img += sino[ib, ia]
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

        img = @. (img + wl * sino[il, ia] + wr * sino[ir, ia])
    end

    img .*= (π * ia_skip / na)
    return embed(img, mask)::Matrix{T}
end
=#


# new version with threads and simplified code

"""
    fbp_back_par(sino, angles,
        ds, offset, xc, yc, mask ; ia_skip, T)
Pixel-driven back-projection
for a grid of `(xc,yc)` pixel center locations
for sinogram `sino` from a
parallel-beam geometry.
It assumes the angles are equally spaced over `[0,π)`.

# in
- `sino::Matrix{<:Number}` `(nb, na)` usually ramp-filtered
- `angles::Vector{<:Real}` `(na)` in radians
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
function fbp_back_par(
    sino::AbstractMatrix{Ts},
    angles::AbstractVector{To},
    ds::Tds,
    offset::Toffset,
    xc::AbstractArray{Tc},
    yc::AbstractArray{Tc},
    mask::AbstractMatrix{Bool} ;
    ia_skip::Int = 1,
    T = typeof(oneunit(Ts) * (oneunit(To) * oneunit(Tc) / oneunit(Tds) + oneunit(Toffset)))
) where {Ts <: Number, To <: RealU, Tds <: RealU, Toffset <: Real, Tc <: RealU}

    image = zeros(T, size(mask)) # need zero(T) outside mask
    sinϕ = sin.(angles[1:ia_skip:end])
    cosϕ = cos.(angles[1:ia_skip:end])
    fbp_back_par!(image, sino, sinϕ, cosϕ,
        ds, offset, xc, yc, mask ; ia_skip)
    return image
end


"""
    fbp_back_par!(image, sino, sinϕ, cosϕ,
        ds, offset, xc, yc, mask ; ia_skip)
Mutating version of
pixel-driven back-projection
for a grid of `(xc,yc)` pixel center locations
for sinogram `sino` from a parallel-beam geometry.
It uses `Threads`.
It assumes the angles are equally spaced over `[0,π)`.

# in
- `sino::Matrix{<:Number}` `(nb, na)` usually ramp-filtered
- `sinϕ::Vector{<:Real}` `(na)`
- `cosϕ::Vector{<:Real}` `(na)`
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
function fbp_back_par!(
    image::AbstractMatrix{T},
    sino::AbstractMatrix{<:Number},
    sinϕ::AbstractVector{<:Real}, # sin.(angles) (na_subset)
    cosϕ::AbstractVector{<:Real}, # cos.(angles) (na_subset)
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
        image[c] = fbp_back_par_xy(
            sino, sinϕ, cosϕ,
            wb, xc_ds[c[1]], yc_ds[c[2]]; T,
        )
    end

#=
    # similar speed
    for iy in 1:size(mask,2), ix in 1:size(mask,1)
        if mask[ix,iy]
            image[ix,iy] = fbp_back_par_xy(
                sino, sinϕ, cosϕ, wb,
                xc_ds[ix], yc_ds[iy]; T,
            )
        end
    end
=#

    return image
end


"""
    fbp_back_par_xy(sino, sinϕ, cosϕ,
        wb, x, y ; T)
Pixel-driven back-projection for a single (x,y) location
for sinogram `sino` from a parallel-beam geometry.
It assumes the angles are equally spaced over `[0,π)`.

# in
- `sino::Matrix{<:Number}` `(nb, na)` usually ramp-filtered
- `sinϕ::Vector{<:Real}` `(na)`
- `cosϕ::Vector{<:Real}` `(na)`
- `wb::Real = (nb+1)/2 + offset` where usually `offset=0`
- `x,y::Real` pixel center location, normalized by ray spacing

# option
- `T::Type{<:Number}` typically same as `eltype(sino)`

# out
- Returns a scalar of type `T`.
"""
function fbp_back_par_xy(
    sino::AbstractMatrix{Ts}, # (nb, na_subset)
    sinϕ::AbstractVector{To}, # sin.(angles) (na_subset)
    cosϕ::AbstractVector{To}, # cos.(angles) (na_subset)
    wb::Tb, # (nb+1)/2 + offset
    x_ds::Tx, # xc / ds
    y_ds::Tx ;
    T::Type{<:Number} = typeof(oneunit(Ts) * one(To) * one(Tb) * one(Tx)),
) where {Ts <: Number, To <: Real, Tb <: Real, Tx <: Real}

    nb = size(sino,1)
    na_subset = size(sino,2)

    pixel = zero(T)

    for ia in 1:na_subset
        @inbounds bb = x_ds * cosϕ[ia] + y_ds * sinϕ[ia] + wb # unitless bin index

#=
        # nearest neighbor interpolation:
        ib = round(Int, bb)
        if 1 ≤ ib ≤ nb
            @inbounds pixel += sino[ib, ia]
        end
=#

        # linear interpolation:
        il = floor(Int, bb) # left bin
        ir = 1 + il # right bin

        if (1 ≤ il) && (ir ≤ nb)
            wr = bb - il # left weight
            wl = 1 - wr # right weight
            @inbounds pixel += wl * sino[il, ia] + wr * sino[ir, ia]
        end
    end

    return pixel # * (π / na_subset) # todo
end

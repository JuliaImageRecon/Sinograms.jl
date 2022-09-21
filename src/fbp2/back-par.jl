# fbp2/back-par.jl

export fbp_back

using LazyGrids: ndgrid
using ImageGeoms: ImageGeom, embed
#using Sinograms: SinoPar



"""
    img = fbp_back(sg, ig, sino ; ia_skip)

2D pixel-driven backprojection for FBP.

in
- `sg::SinoGeom`
- `ig::ImageGeom`
- `sino::AbstractArray{<:Number}` sinogram(s) (line integrals), usually ramp filtered

options
- `ia_skip::Int` downsample in angle to save time for quick tests (default: 1)

out
- `img::Array{<:Number}` reconstructed image(s)
"""
fbp_back

# parallel-beam case
function fbp_back(
    sg::SinoPar{Td, To},
    ig::ImageGeom,
    sino::AbstractMatrix{Ts} ;
    ia_skip::Int = 1,
#   do_r_mask::Bool = false
) where {Td, To, Ts <: Number}

    sg.dim == size(sino) || throw("sino size")

    # type inference help:
    Toffset = Float32 # eltype(sg.offset)
    T = eltype(oneunit(Ts) * (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

    return fbp_back_par(
        sino, sg.ar,
        sg.ds, sg.offset,
        sg.rfov,
        ndgrid(axes(ig)...)...,
        ig.mask, ia_skip,
    )::Matrix{T}
end


function fbp_back_par(
    sino::AbstractMatrix{Ts},
    angles::AbstractVector{To},
    ds::Tds,
    offset::Toffset,
    rfov::RealU,
    xc::AbstractArray{Tc},
    yc::AbstractArray{Tc},
    mask::AbstractMatrix{Bool},
    ia_skip::Int,
) where {Tds <: RealU, Tc <: RealU, Toffset <: Real, Ts <: Number, To <: RealU}

    Td = promote_type(Tds, Tc)
    T = eltype(oneunit(Ts) * (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

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

    warned = false
    for ia in 1:ia_skip:na

        angle = angles[ia]
        (sang, cang) = sincos(angle)
        bb = @. (xc * cang + yc * sang) # [np]
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


    # trick: extra zero column saves linear interpolation indexing within loop!
#   sino_padded = cat(dims=1, sino, zeros(Ts, 1, size(sino,2)))
#   sino_padded[:,1:ia_skip:end]
#   sin.(angles[1:ia_skip:end])
#   cos.(angles[1:ia_skip:end])


#=
# back-projection for a single (x,y) location (todo)
function fbp_back_par_xy(
    sino_padded::AbstractMatrix{Ts}, # (nb+1, na_subset)
    sang::AbstractVector{To}, # sin.(angles) (na_subset)
    cang::AbstractVector{To}, # cos.(angles) (na_subset)
    ds::Td,
    offset::Toffset,
    rfov::RealU,
    x::Td,
    y::Td ;
    warn_truncated::Bool = true,
) where {Td <: RealU, Toffset <: Real, Ts <: Number, To <: RealU}

    T = eltype(oneunit(Ts) * (oneunit(Td) * oneunit(To) / oneunit(Td) + oneunit(Toffset)))

    nb1, na_subset = size(sino_padded)
    nb = nb1 - 1

    wb = Toffset((nb+1)/2 + offset)

    img = zero(T)

    bb = @. (x * cang + y * sang) # [na or a subset thereof]
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
        warn_truncated && @warn("image exceeds system FOV; modify ig.mask?")
        il[.!ig] .= nb+1
        ir[.!ig] .= nb+1
    end

    wr = bb - il # left weight
    wl = 1 .- wr # right weight

    fun(ia) = wl[ia] * sino[il[ia], ia] + wr[ia] * sino[ir[ir], ia]
    return sum(fun, 1:na_subset) * (π / na_subset)
end
=#

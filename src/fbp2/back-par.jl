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
    sg::SinoPar,
    ig::ImageGeom,
    sino::AbstractMatrix{<:Number} ;
    ia_skip::Int = 1,
#   do_r_mask::Bool = false
)

    sg.dim == size(sino) || throw("sino size")

    return fbp_back_par(
        sino, sg.ar,
        sg.ds, sg.offset,
        sg.rfov,
        ndgrid(axes(ig)...)...,
        ig.mask, ia_skip,
    )
end


function fbp_back_par(
    sino::AbstractMatrix{<:Number},
    angles::AbstractVector,
    ds::RealU,
    offset::Real,
    rfov::RealU,
    xc,
    yc,
    mask::AbstractMatrix{Bool},
    ia_skip::Int,
)

    nb, na = size(sino)

    # trick: extra zero column saves linear interpolation indexing within loop!
    sino = cat(dims=1, sino, zeros(eltype(sino), 1, size(sino,2)))

#=
    if do_r_mask
        rr = @. sqrt(abs2(xc) + abs2(yc)) # (nx,ny)
        mask = mask .& (rr .< rfov)
    end
=#

    xc = xc[vec(mask)] # [np] pixels within mask
    yc = yc[vec(mask)]

    wb = (nb+1)/2 + offset

    img = 0

    for ia in 1:ia_skip:na

        angle = angles[ia]
        (sang, cang) = sincos(angle)
        bb = @. (xc * cang + yc * sang) # [np]
        bb = @. (bb / ds + wb + 1) # unitless bin index, +1 because julia

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

        # deal with truncated sinograms
        ig = (il .≥ 1) .& (ir .≤ nb)
        il[.!ig] .= nb+1
        ir[.!ig] .= nb+1

#       if !do_r_mask
#       il = min.(il, nb+1)
#       il = max.(il, 1)
#       end

        wr = bb - il # left weight
        wl = 1 .- wr # right weight

        img = @. (img + wl * sino[il, ia] + wr * sino[ir, ia])
    end

    img .*= (π * ia_skip / na)
    return embed(img, mask)
end

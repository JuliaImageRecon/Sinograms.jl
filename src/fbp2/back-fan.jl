# fbp2/back-fan.jl

export fbp_back

using LazyGrids: ndgrid
using ImageGeoms: ImageGeom, embed
#using Sinograms: SinoFan


# fan-beam case
function fbp_back(
    sg::SinoFan,
    ig::ImageGeom,
    sino::AbstractMatrix{<:Number} ;
    ia_skip::Int = 1,
)

    sg.dim == size(sino) || throw("sino size")

    is_arc = sg.dfs == 0 ? true : isinf(sg.dfs) ? false : throw("bad dfs")

    return fbp_back_fan(
        sino, sg.ar,
        sg.dsd, sg.dso, sg.dfs,
        sg.source_offset, is_arc,
        sg.ds, sg.offset,
        sg.rfov,
        ndgrid(axes(ig)...)...,
        ig.mask, ia_skip,
    )
end


function fbp_back_fan(
    sino::AbstractMatrix{<:Number},
    betas::AbstractVector,
    dsd::RealU,
    dso::RealU,
    dfs::RealU,
    source_offset::RealU,
    is_arc::Bool,
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

        # deal with truncated sinograms
        ig = (il .≥ 1) .& (ir .≤ nb)
        il[.!ig] .= nb+1
        ir[.!ig] .= nb+1
   #    if any(il < 1 | il >= nb), error 'bug', end

        wr = bb - il # left weight
        wl = 1 .- wr # right weight

        img = @. (img + (wl * sino[il, ia] + wr * sino[ir, ia]) * w2)
    end

    img .*= (π * ia_skip / na)
    return embed(img, mask)
end

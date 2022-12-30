# fbp2/fbp.jl

export fbp

using ImageGeoms: ImageGeom
# using Sinograms: SinoGeom, SinoPar, SinoFan, SinoMoj, _reale, FBPplan


"""
    image, sino_filt = fbp(plan, sino)

Filtered back-projection (FBP) reconstruction,
returning image and filtered sinogram.

in
- `plan::FBPplan`
- `sino::AbstractMatrix{<:Number}`

out
- `image::Matrix{<:Number}`       reconstructed image(s)
- `sino_filt::Matrix{<:Number}`   filtered sinogram(s)
"""
fbp

function fbp(
    plan::FBPNormalPlan{<:SinoPar},
    sino::AbstractMatrix{<:Number},
)
    return fbp(plan.rg, plan.ig, sino, plan.filter, plan.view_weight)
end


function fbp(
    plan::FBPNormalPlan{<:SinoFan},
    sino::AbstractMatrix{<:Number},
)
    sino = sino .* fbp_sino_weight(plan.rg) # fan-beam weighting
    out = fbp(plan.rg, plan.ig, sino, plan.filter, plan.view_weight)
    return out[1], out[2]
end


# 3D stack of sinograms
function fbp(
    plan::FBPNormalPlan{<:SinoGeom},
    aa::AbstractArray{Ts, D},
) where {Ts <: Number, D}
    Th = eltype(plan.filter)
    Tp = eltype(plan.view_weight)
    To = typeof(oneunit(Ts) * oneunit(Th) * oneunit(Tp))
    fun(sino2) = fbp(plan, sino2)[1] # discard sino_filt!
    out = mapslices(fun, aa, dims = [1,2])
    return out::Array{To,D} # todo
end


# FBP for geometries that do not need angle-dependent filters (all but Moj)
function fbp(
    rg::SinoGeom,
    ig::ImageGeom,
    sino::AbstractMatrix{<:Number},
    filter::AbstractVector{<:Number},
    view_weight::AbstractMatrix{<:Number} = ones(1,1),
)
    dims(rg) == size(sino) || error("bad sino size")

    sino_filt = sino .* view_weight
    sino_filt = fbp_sino_filter(sino_filt, filter)

    image = fbp_back(rg, ig, sino_filt)
    return image, sino_filt
end


#=
function fbp(
    plan::FBPNormalPlan{<:SinoFan},
    sino::AbstractMatrix{<:Number},
)

    dfs = plan.rg.dfs
    dfs != 0 && ~isinf(dfs) && throw("only arc or flat fan done")
=#

#=
function fbp(
    plan::FBPNormalPlan{<:SinoMoj}, # todo
    sino::AbstractMatrix{<:Number},
)

    sino = fbp2_apply_sino_filter_moj(sino, plan.moj.H)

    if plan.rg.dx == abs(plan.ig.dx)
        image = plan.moj.G' * sino # backproject
        image = image * (pi / plan.rg.na) # account for "dphi" in integral
    else # revert to conventional pixel driven
        ig = plan.ig
        rg = plan.rg
        arg1 = [uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, sign(ig.dy) * ig.offset_y] # trick: old backproject
        arg2 = [rg.d(1:rg.na), rg.offset, rg.orbit, rg.orbit_start]
# todo  image = jf_mex("back2", arg1[:], arg2[:], int32(arg.nthread), single(sino))
        image = image .* plan.ig.mask
    end

    return image, sino
end
=#


# fbp2_recon_moj
#=
function fpb2(plan::MojPlan, sino::AbstractMatrix{<:Number})
    moj = plan.moj
    # bin sinogram to mojette sampling, whether it is fan or parallel
    msino = moj.ob_rebin * sino
    msino = fbp2_apply_sino_filter_moj(msino, moj.H) # filter mojette sinogram
    image = moj.G' * msino # backproject
    return image * (π / size(sino,2)) # account for "dphi" in integral
end
=#

#=
# todo: needs broadcasts?
function fbp_make_sino_filter_moj(nb, na, dx, orbit, orbit_start, window)
    ang = deg2rad(orbit_start .+ (0:na-1)./na .* orbit)
    npad = nextpow(2, rg.nb + 1) # padded size


    dr = dx * max(abs(cos(ang)), abs(sin(ang)))
    if true
        junk, H = fbp2_sino_filter(:flat, ones(nb,1), ds=1, window=:window, decon1=0)

        # trick: ramp filter usually has 1/dr^2 in it, but convolution by sum
        # requires "dr" so we need dr / dr^2 = 1 / dr
        return H * (1 ./ dr) # [npad,na]

    else # trick: make sure cutoff frequencies match
        u0 = 1/2/dx # standard cutoff for coarsest sampling
        r = -(npad/2):(npad/2-1) * dr # [nb na]
        h = u0^2 .* (2 .* sinc.(2*u0*r) - sinc.(u0*r).^2)
        h = h .* repeat(dr, [npad 1]) # extra dr for discrete-space convolution
   #    clf, plot(r, h, '.'), keyboard
        return _reale(fft(fftshift(h,1), 1))
   #    H = fbp_apodize(H, nb, window)
        !isempty(window) && throw("window not done yet due to dr")
    end
end


function fbp2_apply_sino_filter_moj(sino, H)
    nb = size(sino,1)
    npad = nextpow(2, rg.nb + 1) # padded size
    sinopad = [sino; zeros(npad-nb,size(sino,2))] # padded sinogram
    sino = _reale(ifft(fft(sinopad, 1) .* H), 1)
    sino = sino[1:nb,:]
end
=#

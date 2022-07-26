# fbp2.jl

export plan_fbp, fbp, FBPplan

using ImageGeoms: ImageGeom
# using Sinograms: SinoGeom, SinoPar, SinoFan, SinoMoj, parker_weight, _reale

abstract type FBPplan end

#=
struct Moj{Th,Tg}
    H::Th
    G::Tg
end

Moj() = Moj{AbstractMatrix{<:Real},AbstractMatrix{<:Real}}(zeros(Real,1,1),zeros(Real,1,1))
#Moj(a,b) = Moj{AbstractMatrix{<:Real},AbstractMatrix{<:Real}}(zeros(Real,1,1),zeros(Real,1,1)) # todo method overwritten warning
=#


struct NormalPlan{
    S <: SinoGeom,
    I <: ImageGeom,
    H <: AbstractVector{<:RealU},
    P <: AbstractVector{<:Real},
} <: FBPplan
    sg::S
    ig::I
    filter::H # frequency response Hk of apodized ramp filter, length npad
    parker_weight::P
#   moj::Moj # todo

#=
    function NormalPlan(
        sg::S,
        ig::I,
        window::W,
        parker_weight::P,
    ) where {S <: SinoGeom, I <: ImageGeom,
        W <: Window, P <: AbstractMatrix{<:Real}}
        return NormalPlan{S,I,W,P}(sg, ig, window, parker_weight) #, Moj())
    end
=#
end


#=
struct DfPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Window
    parker_weight::AbstractMatrix{<:Real}
end

struct MojPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Window
    parker_weight::AbstractMatrix{<:Real}
    moj::Moj
end

struct TabPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Window
    parker_weight::AbstractMatrix{<:Real}
end
=#


"""
    plan = plan_fbp(sg, ig; how=:normal, window=Window())

Plan FBP 2D tomographic image reconstruction
for parallel-beam & fan-beam cases,
with either flat or arc detector for fan-beam case.

To use this, you first call it with the sinogram and image geometry.
The routine returns the initialized "plan".
Thereafter, to to perform FBP reconstruction,
call `fbp` with the plan
(perhaps numerous times for the same geometry).


# in
- `sg::SinoGeom`
- `ig::ImageGeom` # only reconstruct pixels within `ig.mask`.

# options
- `how::Symbol`  how to reconstruct
    * `:normal`  default
    * `:mojette` use mojette rebinning and Gtomo2_table
- `window::Window` e.g. `Window(Hamming(0.5))`; default `Window()`
= `npad::Int` # of radial bins after padding; default `nextpow(2, sg.nb + 1)`
-`T::DataType`              type of sino elements (default: `Float32`)

# out
- `plan::FBPplan`            initialized plan

"""
function plan_fbp(
    sg::SinoGeom = SinoPar(),
    ig::ImageGeom = ImageGeom() ;
    how::Symbol = :normal,
    window::Window = Window(),
    npad::Int = nextpow(2, sg.nb + 1),
    decon1::Bool = true,
    T::DataType = Float32,
#   nthread::Int = Threads.nthreads()
)

    weight = parker_weight(sg)
    filter = fbp_filter(sg ; npad, window, decon1)

#   if how === :normal
        return NormalPlan(sg, ig, filter, weight)
#       plan = plan_fbp_normal(sg, ig, filter)
#       plan = plan_fbp_normal(sg, ig, npad, window, filter ; T)
#   elseif how === :dsc
#       plan = plan_fbp_dsc(sg, ig, how, window)
#   elseif how === :df
#       plan = plan_fbp_df(sg, ig, how, window)
#   elseif how === :mojette
#       plan = plan_fbp_moj(sg, ig, how, window)
#   elseif how === :table
#       plan = plan_fbp_tab(sg, ig, how, window)
#   else
#       throw("unknown type: $how")
#   end

#   return plan
end


#=
function plan_fbp_normal(
    sg::SinoGeom,
    ig::ImageGeom,
#   npad::Int,
#   window::Window,
    filter::AbstractVector{<:RealU},
#   T::DataType = Float32,
)
    weight = parker_weight(sg)

#   return NormalPlan(sg, ig, window, weight, filter, npad)
    return NormalPlan(sg, ig, weight, filter)
end
=#


#=
function plan_fbp_normal(
    sg::SinoMoj,
    ig::ImageGeom ;
    window::Window = Window(),
    T::DataType = Float32,
)
    weight = ones(T, sg.nb, sg.na)
    if sg.dx == abs(ig.dx)
        throw("not done")
#       plan.moj.G = Gtomo2_table(sg, ig, ["mojette,back1"], nthread=nthread)
    else
        d = sg.dx
        dx = ig.dx
        @warn("mojette sinogram with d=$d vs image with dx=$dx")
    end

    plan.moj.H = fbp2_make_sino_filter_moj(sg.nb, sg.na, sg.dx, sg.orbit, sg.orbit_start, window)
    return NormalPlan{SinoMoj}(sg, ig, window, weight, moj)
end
=#


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
function fbp(
    plan::NormalPlan{<:SinoGeom},
    sino::AbstractMatrix{<:Number},
)
    return fbp(plan.sg, plan.ig, sino, plan.filter, plan.parker_weight)
end

function fbp(
    plan::NormalPlan{<:SinoFan},
    sino::AbstractMatrix{<:Number},
)
    sino = sino .* fbp_sino_weight(plan.sg) # fan-beam weighting
    return fbp(plan.sg, plan.ig, sino, plan.filter, plan.parker_weight)
end


# 3D stack of sinograms
function fbp(
    plan::NormalPlan{<:SinoGeom},
    aa::AbstractArray{<:Number},
)
    fun(sino2) = fbp(plan, sino2)[1] # discard sino_filt!
    return mapslices(fun, aa, dims = [1,2])
end


# FBP for geometries that do not need angle-dependent filters (all but Moj)
function fbp(
    sg::SinoGeom,
    ig::ImageGeom,
    sino::AbstractMatrix{<:Number},
    filter::AbstractVector{<:Number},
    parker_weight::AbstractVector{<:Number} = ones(size(sino,2)),
)
    sg.dim != size(sino) && throw("bad sino size")
    length(parker_weight) == sg.na || throw("bad parker size")

    sino_filt = sino .* parker_weight'
    sino_filt = fbp_sino_filter(sino_filt, filter) # todo: extra?

    image = fbp_back(sg, ig, sino_filt)
    return image, sino_filt
end


#=
function fbp(
    plan::NormalPlan{<:SinoFan},
    sino::AbstractMatrix{<:Number},
)

    dfs = plan.sg.dfs
    dfs != 0 && ~isinf(dfs) && throw("only arc or flat fan done")
=#

#=
function fbp(
    plan::NormalPlan{<:SinoMoj}, # todo
    sino::AbstractMatrix{<:Number},
)

    sino = fbp2_apply_sino_filter_moj(sino, plan.moj.H)

    if plan.sg.dx == abs(plan.ig.dx)
        image = plan.moj.G' * sino # backproject
        image = image * (pi / plan.sg.na) # account for "dphi" in integral
    else # revert to conventional pixel driven
        ig = plan.ig
        sg = plan.sg
        arg1 = [uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, sign(ig.dy) * ig.offset_y] # trick: old backproject
        arg2 = [sg.d(1:sg.na), sg.offset, sg.orbit, sg.orbit_start]
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
    npad = 2^ceil(log2(2*nb-1)) # padded size

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
    npad = 2^ceil(log2(2*nb-1)) # padded size
    sinopad = [sino; zeros(npad-nb,size(sino,2))] # padded sinogram
    sino = _reale(ifft(fft(sinopad, 1) .* H), 1)
    sino = sino[1:nb,:]
end
=#

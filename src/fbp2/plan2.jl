# fbp2/plan2.jl

export FBPplan, plan_fbp

using ImageGeoms: ImageGeom
# using Sinograms: SinoGeom, SinoPar, SinoFan, SinoMoj, parker_weight

"""
    FBPplan
Abstract type for FBP plans.
"""
abstract type FBPplan end

#=
struct Moj{Th,Tg}
    H::Th
    G::Tg
end

Moj() = Moj{AbstractMatrix{<:Real},AbstractMatrix{<:Real}}(zeros(Real,1,1),zeros(Real,1,1))
#Moj(a,b) = Moj{AbstractMatrix{<:Real},AbstractMatrix{<:Real}}(zeros(Real,1,1),zeros(Real,1,1)) # todo method overwritten warning
=#


"""
    FBPNormalPlan{S,I,H,P}
Struct type for storing "normal" FBP plan.
"""
struct FBPNormalPlan{
    S <: SinoGeom,
    I <: ImageGeom,
    H <: AbstractVector{<:RealU},
    P <: Any,
} <: FBPplan
    sg::S
    ig::I
    filter::H # frequency response Hk of apodized ramp filter, length npad
    parker_weight::P # typically a Matrix of nonnegative reals
#   moj::Moj # todo

#=
    function FBPNormalPlan(
        sg::S,
        ig::I,
        window::W,
        parker_weight::P,
    ) where {S <: SinoGeom, I <: ImageGeom,
        W <: Window, P <: AbstractMatrix{<:Real}}
        return FBPNormalPlan{S,I,W,P}(sg, ig, window, parker_weight) #, Moj())
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
The routine returns the initialized `plan`.
Thereafter, to to perform FBP reconstruction,
call `fbp` with the `plan`
(perhaps numerous times for the same geometry).


# in
- `sg::SinoGeom`
- `ig::ImageGeom` only reconstruct pixels within `ig.mask`.

# options
- `how::Symbol` how to reconstruct
    * `:normal` default
    * `:mojette` use mojette rebinning and Gtomo2_table
- `window::Window` e.g., `Window(Hamming(), 0.8)`; default `Window()`
- `npad::Int` # of radial bins after padding; default `nextpow(2, sg.nb + 1)`
- `decon1::Bool` deconvolve interpolator effect? (default `true`)
- `T::Type{<:Number}` type of `sino` elements (default `Float32`)

# out
- `plan::FBPplan` initialized plan

"""
function plan_fbp(
    sg::SinoGeom,
    ig::ImageGeom ;
    how::Symbol = :normal,
    window::Window = Window(),
    npad::Int = nextpow(2, sg.nb + 1),
    decon1::Bool = true,
    T::Type{<:Number} = Float32,
#   nthread::Int = Threads.nthreads(),
)

    weight = parker_weight(sg)
    filter = fbp_filter(sg ; npad, window, decon1)

#   if how === :normal
        return FBPNormalPlan(sg, ig, filter, weight)
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
#   T::Type{<:Number} = Float32,
)
    weight = parker_weight(sg)

#   return FBPNormalPlan(sg, ig, window, weight, filter, npad)
    return FBPNormalPlan(sg, ig, weight, filter)
end
=#


#=
function plan_fbp_normal(
    sg::SinoMoj,
    ig::ImageGeom ;
    window::Window = Window(),
    T::Type{<:Number} = Float32,
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
    return FBPNormalPlan{SinoMoj}(sg, ig, window, weight, moj)
end
=#


function Base.show(io::IO, ::MIME"text/plain", p::FBPNormalPlan{S,I,H,P}) where {S,I,H,P}
    println(io, "FBPNormalPlan{S,I,H,P} with")
    println(io, " S = $S")
    println(io, " I = $I")
    println(io, " H = $H with extrema ", extrema(p.filter))
    println(io, " P = $P ", size(p.parker_weight), " with extrema ", extrema(p.parker_weight))
end

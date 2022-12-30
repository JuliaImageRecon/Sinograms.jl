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

The `view_weight` can include (products of)
- Parker weighting for short scans
- view-wise weighting from `fbp_sino_weight` todo
- `dβ` weighting for possibly nonuniform angles from `_angle_weights`
"""
struct FBPNormalPlan{
    S <: SinoGeom,
    I <: ImageGeom,
    H <: AbstractVector{<:RealU},
#   P <: Any,
    V <: AbstractArray{<:RealU},
} <: FBPplan
    rg::S
    ig::I
    filter::H # frequency response Hk of apodized ramp filter, length npad
#   parker_weight::P # typically a Matrix of nonnegative reals
    view_weight::V
#   moj::Moj # todo
end


#=
struct DfPlan <: FBPplan
    rg::SinoGeom
    ig::ImageGeom
    window::Window
    parker_weight::AbstractMatrix{<:Real}
end

struct MojPlan <: FBPplan
    rg::SinoGeom
    ig::ImageGeom
    window::Window
    parker_weight::AbstractMatrix{<:Real}
    moj::Moj
end

struct TabPlan <: FBPplan
    rg::SinoGeom
    ig::ImageGeom
    window::Window
    parker_weight::AbstractMatrix{<:Real}
end
=#


"""
    plan = plan_fbp(rg, ig; how=:normal, window=Window())

Plan FBP 2D tomographic image reconstruction
for parallel-beam & fan-beam cases,
with either flat or arc detector for fan-beam case.

To use this method,
you first call it with the sinogram geometry
and image geometry.
The routine returns the initialized `plan`.
Thereafter, to to perform FBP reconstruction,
call `fbp` with the `plan`
(perhaps numerous times for the same geometry).

# in
- `rg::SinoGeom`
- `ig::ImageGeom` only reconstruct pixels within `ig.mask`.

# options
- `how::Symbol` how to reconstruct
    * `:normal` default
    * `:mojette` use mojette rebinning and Gtomo2_table
- `window::Window` e.g., `Window(Hamming(), 0.8)`; default `Window()`
- `npad::Int` # of radial bins after padding;
  default `nextpow(2, rg.nb + 1)`
- `decon1::Bool` deconvolve interpolator effect? (default `true`)
- `T::Type{<:Number}` type of `sino` elements (default `Float32`)

# out
- `plan::FBPplan` initialized plan

"""
function plan_fbp(
    rg::SinoGeom,
    ig::ImageGeom,
    ;
    how::Symbol = :normal,
    window::Window = Window(),
    npad::Int = nextpow(2, rg.nb + 1),
    decon1::Bool = true,
    filter::AbstractVector{<:RealU} = fbp_filter(rg ; npad, window, decon1),
    T::Type{<:Number} = Float32,
#   nthread::Int = Threads.nthreads(),
)

    weight = _fbp_weights(rg)

#   if how === :normal
        return FBPNormalPlan(rg, ig, filter, weight)
#       plan = plan_fbp_normal(rg, ig, filter)
#       plan = plan_fbp_normal(rg, ig, npad, window, filter ; T)
#   elseif how === :dsc
#       plan = plan_fbp_dsc(rg, ig, how, window)
#   elseif how === :df
#       plan = plan_fbp_df(rg, ig, how, window)
#   elseif how === :mojette
#       plan = plan_fbp_moj(rg, ig, how, window)
#   elseif how === :table
#       plan = plan_fbp_tab(rg, ig, how, window)
#   else
#       throw("unknown type: $how")
#   end

#   return plan
end


function _fbp_weights(rg::SinoGeom)
    weight = parker_weight(rg) .*
        _angle_weights(_ar(rg))
    return weight
end


#=
function plan_fbp_normal(
    rg::SinoGeom,
    ig::ImageGeom,
#   npad::Int,
#   window::Window,
    filter::AbstractVector{<:RealU},
#   T::Type{<:Number} = Float32,
)
    weight = parker_weight(rg)

#   return FBPNormalPlan(rg, ig, window, weight, filter, npad)
    return FBPNormalPlan(rg, ig, weight, filter)
end
=#


#=
function plan_fbp_normal(
    rg::SinoMoj,
    ig::ImageGeom ;
    window::Window = Window(),
    T::Type{<:Number} = Float32,
)
    weight = ones(T, rg.nb, rg.na)
    if rg.d == abs(ig.dx)
        throw("not done")
#       plan.moj.G = Gtomo2_table(rg, ig, ["mojette,back1"], nthread=nthread)
    else
        d = rg.d
        dx = ig.dx
        @warn("mojette sinogram with d=$d vs image with dx=$dx")
    end

    plan.moj.H = fbp2_make_sino_filter_moj(rg.nb, rg.na, rg.d, rg.orbit, rg.orbit_start, window)
    return FBPNormalPlan{SinoMoj}(rg, ig, window, weight, moj)
end
=#


function Base.show(io::IO, ::MIME"text/plain", p::FBPNormalPlan{S,I,H,V}) where {S,I,H,V}
    rg = p.rg
    println(io, "FBPNormalPlan{S,I,H,V} with")
    println(io, " S = $S ", (rg.nb, rg.na))
    ig = p.ig
    println(io, " I = $I ", ig.dims)
    println(io, " H = $H with extrema ", extrema(p.filter))
#   println(io, " P = $P ", size(p.parker_weight), " with extrema ", extrema(p.parker_weight))
    println(io, " V = $V ", size(p.view_weight), " with extrema ", extrema(p.view_weight))
end

# fbp3/plan3.jl

export FDKplan, plan_fbp

using ImageGeoms: ImageGeom
# using Sinograms: CtGeom Window fbp_filter parker_weight_fan


"""
    FDKplan{C,I,H,P}
Struct type for storing FDK plan.
"""
struct FDKplan{
    C <: CtGeom,
    I <: ImageGeom{3},
    H <: AbstractVector{<:RealU},
    P <: Any,
}
    cg::C
    ig::I
    filter::H # frequency response Hk of apodized ramp filter, length npad
    parker_weight::P

#=
    function FDKplan(
        cg::C,
        ig::I,
        window::W,
        parker_weight::P,
    ) where {C <: CtGeom, I <: ImageGeom{3},
        W <: Window, P <: AbstractArray{<:Real}}
        return FDKplan{C,I,W,P}(cg, ig, window, parker_weight)
    end
=#
end


function parker_weight(cg::CtFan; kwargs...)
    return parker_weight_fan(
        cg.ns, cg.na, cg.orbit, cg.orbit_short,
        cg.ar, cg.gamma, cg.gamma_max; kwargs...,
    )
end


"""
    plan = plan_fbp(cg, ig; window=Window(), ...)

Plan FDK 3D CBCT image reconstruction,
with either flat or arc detector.

To use this, you first call it with the CT geometry and image geometry.
The routine returns the initialized `plan`.
Thereafter, to to perform FDK reconstruction,
call `fbp` with the `plan`
(perhaps numerous times for the same geometry).


# in
- `cg::CtGeom`
- `ig::ImageGeom` only reconstruct pixels within `ig.mask`.

# options
- `window::Window` e.g., `Window(Hamming(), 0.8)`; default `Window()`
- `npad::Int` # of radial bins after padding; default `nextpow(2, cg.ns + 1)`
- `decon1::Bool` deconvolve interpolator effect? (default `true`)
-`T::Type{<:Number}` type of sino elements (default `Float32`)

# out
- `plan::FDKplan` initialized plan

"""
function plan_fbp(
    cg::CtGeom,
    ig::ImageGeom ;
    window::Window = Window(),
    npad::Int = nextpow(2, cg.ns + 1),
    decon1::Bool = true,
    T::Type{<:Number} = Float32,
)

    filter = fbp_filter(cg ; npad, window, decon1)
    weight = parker_weight(cg)

    return FDKplan(cg, ig, filter, weight)
end


function Base.show(io::IO, ::MIME"text/plain", p::FDKplan{C,I,H,P}) where {C,I,H,P}
    c = p.cg
    println(io, "FDKplan{C,I,H,P} with")
    println(io, " C = $C ", (c.ns, c.nt, c.na))
    i = p.ig
    println(io, " I = $I ", i.dims)
    println(io, " H = $H ", size(p.filter), " with extrema ", extrema(p.filter))
    println(io, " P = $P ", size(p.parker_weight), " with extrema ", extrema(p.parker_weight))
end

# fbp3/plan3.jl

export FDKplan, plan_fbp

using ImageGeoms: ImageGeom
# using Sinograms: CtGeom Window fbp_filter
# using Sinograms: parker_weight _view_weights fdk_weight_cyl


"""
    FDKplan{C,I,H,V}
Struct type for storing FDK plan.

The `view_weight` can include (products of)
- Parker weighting for short scans
- view-wise CBCT weighting from `fdk_weight_cyl`
- `dβ` weighting for possibly nonuniform angles
"""
struct FDKplan{
    C <: CtGeom,
    I <: ImageGeom{3},
    H <: AbstractVector{<:RealU},
    V <: AbstractArray{<:RealU},
}
    rg::C
    ig::I
    filter::H # frequency response Hk of apodized ramp filter, length npad
    view_weight::V
end


"""
    plan = plan_fbp(rg, ig; window=Window(), ...)

Plan FDK 3D CBCT image reconstruction,
with either flat or arc detector.

To use this method,
you first call it with the CT geometry and image geometry.
The routine returns the initialized `plan`.
Thereafter, to to perform FDK reconstruction,
call `fbp` with the `plan`
(perhaps numerous times for the same geometry).

# in
- `rg::CtGeom`
- `ig::ImageGeom` only reconstruct pixels within `ig.mask`.

# options
- `window::Window` e.g., `Window(Hamming(), 0.8)`; default `Window()`
- `npad::Int` # of radial bins after padding; default `nextpow(2, rg.ns + 1)`
- `decon1::Bool` deconvolve interpolator effect? (default `true`)

# out
- `plan::FDKplan` initialized plan

"""
function plan_fbp(
    rg::CtGeom,
    ig::ImageGeom,
    ;
    window::Window = Window(),
    npad::Int = nextpow(2, rg.ns + 1),
    decon1::Bool = true,
    filter::AbstractVector{<:RealU} = fbp_filter(rg ; npad, window, decon1),
)

    weight = _fdk_weights(rg)
    return FDKplan(rg, ig, filter, weight)
end


function _fdk_weights(rg::CtGeom)
    weight = parker_weight(rg)
    weight = weight .* _view_weights(_ar(rg))
    weight = weight .* fdk_weight_cyl(rg)
    return weight
end


function Base.show(io::IO, ::MIME"text/plain", p::FDKplan{C,I,H,V}) where {C,I,H,V}
    rg = p.rg
    println(io, "FDKplan{C,I,H,V} with")
    println(io, " C = $C ", (rg.ns, rg.nt, rg.na))
    i = p.ig
    println(io, " I = $I ", i.dims)
    println(io, " H = $H ", size(p.filter), " with extrema ", extrema(p.filter))
    println(io, " V = $V ", size(p.view_weight), " with extrema ", extrema(p.view_weight))
end

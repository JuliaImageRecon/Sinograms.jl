#=
fbp3/fdk.jl
Translated from feldkamp.m in MIRT
Copyright 2022-5-11 Jason Hu and Jeff Fessler, University of Michigan
=#

export fdk


# helpers

"""
    fdk_weight_cyl
FDK projection weighting providing "exact" CBCT reconstruction
for cylindrical-like objects that satisfy
``f(x,y,z) = f(x,y,0) ∀z``.
The output is a `(ns,nt)` matrix.
"""
fdk_weight_cyl

function fdk_weight_cyl_arc(
    s::RealU,
    t::RealU,
    dsd::RealU,
    dso::RealU,
)
    return (dso/dsd) * cos(s / dsd) / sqrt(abs2(t / dsd) + 1)
end

function fdk_weight_cyl_flat(
    s::RealU,
    t::RealU,
    dsd::RealU,
    dso::RealU,
)
    return dso / sqrt(abs2(s) + abs2(t) + abs2(dsd))
end

fdk_weight_cyl(rg::CtFanArc) =
    fdk_weight_cyl_arc.(_s(rg), _t(rg)', rg.dsd, _dso(rg))
#   Iterators.map((rg) -> fdk_weight_cyl_arc(rg..., rg.dsd, rg.dso),
#       Iterators.product(ct_geom_s(rg), ct_geom_t(rg)))

fdk_weight_cyl(rg::CtFanFlat) =
    fdk_weight_cyl_flat.(_s(rg), _t(rg)', rg.dsd, _dso(rg))


"""
    image = fdk(plan, proj)

Reconstruct 3D image
from cone-beam computed tomography (CBCT) data
collected with a circular source trajectory
via FDK method.

# in
* `plan::FDKplan`
* `proj` (ns,nt,na) projection views

# out
* `image` (nx,ny,nz) reconstructed image

References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
"""
function fdk(plan::FDKplan, proj::AbstractArray{<:Number,3})
    rg = plan.rg

    size(proj) == dims(plan.rg) ||
        error("size mismatch $(size(proj)) $(dims(plan.rg))")

    # step 1: apply cone-beam weights
    proj = proj .* plan.view_weight # todo: before or after filtering?

    # step 2: filter each projection view
    proj = fbp_sino_filter(proj, plan.filter)

    # step 3: cone-beam backprojection
    image = cbct_back(proj, rg, plan.ig)

    return image
end

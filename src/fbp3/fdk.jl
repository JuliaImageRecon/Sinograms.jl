#=
fbp3/fdk.jl
Translated from feldkamp.m in MIRT
Copyright 2022-5-11 Jason Hu and Jeff Fessler, University of Michigan
=#

export fdk


# helpers

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

# (ns,nt) matrix
fdk_weight_cyl(cg::CtFanArc) =
    fdk_weight_cyl_arc.(cg.s, cg.t', cg.dsd, cg.dso)

fdk_weight_cyl(cg::CtFanFlat) =
    fdk_weight_cyl_flat.(cg.s, cg.t', cg.dsd, cg.dso)


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
    cg = plan.cg

    # step 1: apply weights
    proj = proj .* fdk_weight_cyl(cg) # todo: precompute with plan!

    # step 2: filter each projection view
    proj = fbp_sino_filter(proj, plan.filter)

    # step 3: cone-beam backprojection
    image = cbct_back(proj, cg, plan.ig)

    return image
end

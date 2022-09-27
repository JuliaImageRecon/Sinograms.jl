#=
fbp3/fdk.jl
Translated from feldkamp.m in MIRT
Copyright 2022-5-11 Jason Hu and Jeff Fessler, University of Michigan
=#

export fdk


#=
# the following two structs are necessary when only testing with the mat files
# todo: relocate

export ct_geom2
export im_geom2
struct ct_geom2
    type
    ns
    nt
    na
    down
    nframe
    frame
    orbit_start
    pitch
    source_z0
    units
    user_source_zs
    orbit
    ds
    dt
    offset_s
    offset_t
    dsd
    dso
    dod
    dfs
end

struct im_geom2
    nx
    ny
    nz
    dx
    dy
    dz
    offset_x
    offset_y
    offset_z
    offsets
    fov
    zfov
    down
    mask
    is3
    dim
end
=#

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


#=
function feldkamp_do(proj, cg, ig, ds, dt, offset_s, offset_t, offset_source, dsd, dso, dfs, orbit, orbit_start, mask, nz, dx, dy, dz, offset_xyz, w1cyl, window, ia_skip, extrapolate_t, use_mex, nthread)
    #step 1 compute the weights
    proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl)

    #step 2 filter each projection view
    (ns, nt, na) = size(proj)
    proj = fbp_filter(proj, window, dsd, dfs, ds)

    #step 3 cone beam backprojection
    img = cbct_back(proj, cg, ig)

    return img, proj
end
=#


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


#=
"""
    image = feldkamp(cg, ig, proj)

FBP reconstruction of cone-beam computed tomography (CBCT) data
collected with a circular source trajectory.
See feldkamp_example.jl for an example.

# in
* `cg::CtGeom`
* `ig::ImageGeom`
* `proj` (ns,nt,na) projection views

# out
* `img` (nx,ny,nz)    reconstructed image
* `proj_out [ns nt na]    filtered projections (for debugging)

References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
"""

function feldkamp(cg::CtGeom, ig::ImageGeom, proj::AbstractArray{3,<:Number})
#=
    if isa(cg, ct_geom2)
        #using the mat file code, so all the variables are already available
        return feldkamp_do(proj, cg, ig, cg.ds, cg.dt, cg.offset_s, cg.offset_t, 0, cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ig.mask, ig.nz, ig.dx, ig.dy, ig.dz, [ig.offset_x, ig.offset_y, ig.offset_z], 1, "ramp", 1, 0, 0, 0)
    elseif isa(cg, CtFanArc)
        println(size(ig.mask_or))
        return feldkamp_do(proj, cg, ig, cg.ds, cg.dt, cg.offset_s, cg.offset_t, 0, cg.dsd, cg.dso, Inf, cg.orbit, cg.orbit_start, ig.mask, ig.nz, ig.dx, ig.dy, ig.dz, [ig.offset_x, ig.offset_y, ig.offset_z], 1, "ramp", 1, 0, 0, 0)
    elseif isa(cg, CtFanFlat)
        return feldkamp_do(proj, cg, ig, cg.ds, cg.dt, cg.offset_s, cg.offset_t, 0, cg.dsd, cg.dso, 0, cg.orbit, cg.orbit_start, ig.mask, ig.nz, ig.dx, ig.dy, ig.dz, [ig.offset_x, ig.offset_y, ig.offset_z], 1, "ramp", 1, 0, 0, 0)
    else
        error("Feldkamp not implemented for this case")
    end
    return 0
=#

    return image
end
=#

#=
fbp3/cb_arc_to_par.jl
Translated from ir_coord_cb_arc_to_par.m in MIRT
2022-05-12, Jason Hu, Jeff Fessler, University of Michigan
=#

"""
    u, v, azim, polar = cb_arc_to_par(s, t, β, dso, dod)

Convert from cone-beam arc (3rd gen CT) coordinates
to parallel-beam coordinates.

# in
* `s,t` detector coordinates
* `β` X-ray source angle, measured counter-clockwise from y axis
* `dso,dsd` distances for the geometry

# out
* `u,v` transaxial and axial parallel-beam detector coordinates
* `azim` ϕ transaxial or azimuthal angle (radians)
* `polar` θ polar angle (radians)
"""
function cb_arc_to_par(
    s::Ts, t::RealU, β::RealU,
    dso::RealU, dod::RealU,
) where {Ts <: RealU}

    dsd = dso + dod
    dfs = zero(Ts)
    Rf = dfs + dsd
    pos_det = (Rf * sin(s / Rf), dso + dfs - Rf * cos(s / Rf), t)
    pos_src = (zero(Ts), dso, zero(Ts)) # todo pitch / source_zs here
    e = pos_det .- pos_src
    e = e ./ sqrt(sum(abs2, e))

    azim = -atan(e[1], e[2])
    polar = asin(e[3]) # note sign!
    sϕ, cϕ = sincos(azim)
    u = cϕ * pos_src[1] + sϕ * pos_src[2]
    v = (sϕ * pos_src[1] - cϕ * pos_src[2]) * sin(polar) +
        pos_src[3] * cos(polar)
    return u, v, azim + β, polar
end

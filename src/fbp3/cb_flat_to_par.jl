#=
fbp3/cb_flat_to_par.jl
Translated from ir_coord_cb_flat_to_par.m in MIRT
2022-05-12, Jason Hu, Jeff Fessler, University of Michigan
=#

"""
    u, v, azim, polar = cb_flat_to_par(s, t, β, dso, dod)

Convert from cone-beam flat panel coordinages
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
function cb_flat_to_par(s::RealU, t::RealU, β::RealU, dso::RealU, dod::RealU)
    dsd = dso + dod
    u = dso * s / sqrt(s^2 .+ dsd^2)
    v = dso * t / sqrt(s^2 + t^2 + dsd^2) * dsd / sqrt(s^2 + dsd^2)
    azim = atan(s / dsd)
    polar = atan(t / sqrt(s^2 .+ dsd^2)) # note sign!
    return u, v, azim + β, polar
end

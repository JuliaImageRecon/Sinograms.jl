#=
fbp3/cb_flat_to_par.jl
Translated from ir_coord_cb_flat_to_par.m in MIRT
2022-05-12, Jason Hu, Jeff Fessler, University of Michigan
=#

"""
    u, v, azim, polar = cb_flat_to_par(s, t, dso, dod)

Convert from cone-beam flat panel to parallel-beam coordinates.

# in
* `s,t` detector coordinates
* `dso,dsd` distances for the geometry

# out
* `u,v` transaxial and axial parallel-beam detector coordinates
* `azim` transaxial or azimuthal angle (radians)
* `polar` polar angle (radians)
"""
function cb_flat_to_par(s::RealU, t::RealU, dso::RealU, dod::RealU)
    dsd = dso + dod
    u = dso * s / sqrt(s^2 .+ dsd^2)
    v = dso * t / sqrt(s^2 + t^2 + dsd^2) * dsd / sqrt(s^2 + dsd^2)

    azim = atan(s / dsd)
    polar = -atan(t / sqrt(s^2 .+ dsd^2))
    return u, v, azim, polar
end

function cb_flat_to_par(s::RealU, t::RealU, β::RealU, dso::RealU, dod::RealU)
   (u, v, ϕ, θ) = cb_flat_to_par(s, t, dso, dod) # for β=0
   ϕ += β # account for β
   (sϕ, cϕ) = sincos(ϕ)
   (u, v) = u * cϕ + v * sϕ, u * sϕ - v * cϕ
   return (u, v, ϕ, θ)
end

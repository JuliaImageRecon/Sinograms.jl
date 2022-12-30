#=
fbp-par.jl
Simple interfaces to parallel-beam FBP for user convenience
=#

using ImageGeoms: ImageGeom, circle
#using Sinograms: plan_fbp, fbp, RealU, SinoPar

export fbp, fbp!


"""
    fbp(sino::AbstractMatrix ; kwargs...)

FBP reconstruction from a parallel-beam sinogram of size `[nr × nϕ]`.
Returns an image of size `[nx × ny]`.

# Options
* `nx` : default `nr`
* `ny` : default `nx`
* `kwargs` : passed to `fbp!`
"""
function fbp(
    sino::AbstractMatrix{Ts} ;
    nx::Int = size(sino, 1),
    ny::Int = nx,
    dr::RealU = 1,
    dx::RealU = dr,
    dy::RealU = dx,
    deltas::NTuple{2,Td} = (dx, dy),
    kwargs...
) where {Ts <: Number, Td <: RealU}
    U = typeof(1f0 * oneunit(Ts) / oneunit(Td))
    image = zeros(U, nx, ny)
    fbp!(image, sino; dr, dx, dy, deltas, kwargs...)
end


"""
    fbp!(image, sino ; orbit::Real = 180, kwargs...)

FBP reconstruction from a parallel-beam sinogram of size `[nr × nϕ]`.
Writes result into `image` matrix.

# Input
* `sino::AbstractMatrix`

# Options for `SinoPar` constructor
* `dr` : sinogram radial spacing; default 1
* `orbit` : angular range in degrees; default 180
* `orbit_start` : angular range in degrees; default 0

# Options for `ImageGeom`
* `dx`, `dy`, `deltas`, `offset_x`, `offset_y`, `offsets`
* `rmax` maximum radius for mask

# Options
* `kwargs` : passed to `plan_fbp`

# Output
* `image::AbstractMatrix` is mutated
"""
function fbp!(
    image::AbstractMatrix{<:Number},
    sino::AbstractMatrix{<:Number} ;
    dr::RealU = 1,
    orbit::RealU = 180,
    orbit_start::RealU = zero(orbit),
    dx::RealU = dr,
    dy::RealU = dx,
    deltas = (dx, dy),
    offset_x::Real = 0,
    offset_y::Real = 0,
    offsets = (offset_x, offset_y),
    rmax::RealU = zero(dr),
    kwargs...
)
    nb, na = size(sino)
    nx, ny = size(image)
    rg = SinoPar( ; nb, na, d = dr, orbit, orbit_start)
    ig = ImageGeom(; dims=(nx, ny), deltas, offsets)
    mask = circle(ig ; r = (rmax == zero(dr)) ? _rfov(rg) : rmax)
    ig = ImageGeom(ig.dims, ig.deltas, ig.offsets, mask)
    plan = plan_fbp(rg, ig ; kwargs...)
#   fbp!(image, sino, plan) # todo
    tmp, _ = fbp(plan, sino) # todo: integrate better!
    copyto!(image, tmp)
    return image
end


#=
function fbp!(
    image::AbstractMatrix{<:Number},
    sino::AbstractMatrix{<:Number},
    plan::FBPNormalPlan,
)
    tmp, _ = fbp(plan, sino) # todo: integrate better!
    copyto!(image, tmp)
    return image
end
=#

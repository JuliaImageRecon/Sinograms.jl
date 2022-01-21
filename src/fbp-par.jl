#=
fbp-par.jl
Simple interfaces to parallel-beam FBP for user convenience
=#

using MIRT: sino_geom, image_geom # todo
#using Sinograms: fbp2, RealU

export fbp, fbp!


"""
    fbp(sino::AbstractMatrix ; kwargs...)

FBP reconstruction from a parallel-beam sinogram of size `[nr × nϕ]`.
Returns an image of size `[nx × ny]`.

# Options
* `nx` : default `nr`
* `ny` : default `nx`
* `kwargs` : passed to `fbp2`
"""
function fbp(
    sino::AbstractMatrix{T} ;
    nx::Int = size(sino, 1),
    ny::Int = nx,
    kwargs...
) where {T <: Number}
    image = zeros(promote_type(T, Float32), nx, ny)
    fbp!(image, sino; kwargs...)
end



"""
    fbp!(image, sino ; orbit::Real = 180, kwargs...)

FBP reconstruction from a parallel-beam sinogram of size `[nr × nϕ]`.
Writes result into `image` matrix.

# Input
* `sino::AbstractMatrix`

# Options for `sino_geom`
* `dr` : sinogram radial spacing; default 1
* `orbit` : angular range in degrees; default 180
* `orbit_start` : angular range in degrees; default 180

# Options for `image_geom`
* `dx`, `dy`, `offset_x`, `offset_y`

# Options
* `kwargs` : passed to `fbp2`

# Output
* `image::AbstractMatrix` is mutated
"""
function fbp!(
    image::AbstractMatrix{<:Number},
    sino::AbstractMatrix{<:Number} ;
    dr::RealU = 1,
    orbit::Real = 180,
    orbit_start::Real = 0,
    dx::RealU = dr,
    dy::RealU = dx,
    offset_x::Real = 0,
    offset_y::Real = 0,
    kwargs...
)
    nb, na = size(sino)
    nx, ny = size(image)
    sg = sino_geom(:par ; nb, na, d = dr, orbit, orbit_start)
    ig = image_geom(; nx, ny, dx, dy, offset_x, offset_y)
#   plan = FBPplan(sg, ig)
    plan = fbp2(sg, ig) # todo FBPplan ??
#   fbp2!(image, sino, plan) # todo
    tmp, _ = fbp2(plan, sino)
    copyto!(image, tmp)
    return image
end

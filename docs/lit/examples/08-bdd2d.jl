#=
# [2D Branchless Distance-driven Projection and Backprojection](@id 08-bdd2d)

This page describes the 2D branchless distance-driven projection
and backprojection (bdd_2d) for fan-beam geometries with a "flat" detector
using the Julia package
[`Sinograms.jl`](https://github.com/JuliaImageRecon/Sinograms.jl).

This page compares the results to the radon transform method. 

This page was generated from a single Julia file:
[08-bdd2d.jl](@__REPO_ROOT_URL__/08-bdd2d.jl).
=#

# ### Setup

# Packages needed here.

using Sinograms: projection, backprojection, SinoFanFlat, rays, plan_fbp, fbp, Window, Hamming, sino_geom_plot!
using ImagePhantoms: shepp_logan, SheppLoganToft, radon, phantom
using ImageGeoms: ImageGeom, fovs, MaskCircle
using MIRTjim: jim, prompt
using Unitful: mm, @u_str

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

#=
### Fan-beam sinogram of SheppLoganToft phantom

For illustration,
we start by synthesizing
a fan-beam sinogram
of the SheppLoganToft phantom.

For completeness,
we use units (from Unitful),
but units are optional.

=#

# Make a tuple to define image geometry.
deg = 1
geo = (DSD = 949mm, DS0 = 541mm, pSize = 1mm, dSize = 1.0239mm, 
    nPix = 512, nDet = 888, angle = deg2rad.(0:deg:360-deg))

# Use `ImageGeom` to define the image geometry for radon test.
ig = ImageGeom(MaskCircle(); dims=(512,512), deltas = (1mm,1mm) )

# Use `SinoFanFlat` to define the sinogram geometry.
rg = SinoFanFlat( ; nb = 888, d = 1.0239mm, na = 360, dsd = 949mm, dod = 408mm)

# Examine the geometry to verify the FOV:
jim(axes(ig), ig.mask; prompt=false)
sino_geom_plot!(rg, ig)

#
prompt()

# Ellipse parameters for SheppLoganToft phantom:
μ = 0.01 / mm # typical linear attenuation coefficient
ob = shepp_logan(SheppLoganToft(); fovs = fovs(ig), u = (1, 1, μ))

# Create a phantom image:
testimage = phantom(axes(ig)..., ob)

# Show the true phantom image.
jim(reverse(testimage, dims=2))

# Arc fan-beam sinogram for SheppLoganToft phantom:
sinogramB = projection(reverse(rot180(testimage'), dims=2), geo)
sinogramR = radon(rays(rg), ob)
p1 = jim(1:geo.nDet, 1:length(geo.angle), sinogramB'; title="bdd_2d sinogram", aspect_ratio = :auto, xlabel="r (mm)", ylabel="ϕ")
p2 = jim(axes(rg), sinogramR; title="Radon test sinogram", xlabel="r", ylabel="ϕ")
jim(p1,p2)

# Here is the difference image of sinogram using bdd_2d versus radon transform method.
jim(axes(rg), sinogramR - sinogramB', title="Error image", aspect_ratio = :auto, xlabel="r", ylabel="ϕ")

#=
## Backprojection of bdd_2d
=#

# Note that backprojected image is not a reconstructed image. Refer to the next section for image reconstruction via FBP.
sinogram = sinogramB*u"mm^-1" # include units
imageB = backprojection(sinogram, geo)
jim(1:geo.nPix, 1:geo.nPix, imageB'; title="bdd_2d backprojection", aspect_ratio = :auto, xlabel="mm", ylabel="mm")

#=
## Image reconstruction via FBP
This section compares the reconstructed image using the bdd_2d sinogram and the radon test sinogram.

Here we start with a "plan",
which would save work if we were reconstructing many images.
For illustration we include `Hamming` window.
=#

plan = plan_fbp(rg, ig; window = Window(Hamming(), 1.0))
fbp_image_b = fbp(plan, sinogramB')
fbp_image_r = fbp(plan, sinogramR)

# A narrow color window is needed to see the soft tissue structures:
clim = (0.04, 1.1) .* μ
r1 = jim(axes(ig), fbp_image_b, "FBP image with bdd_2d"; clim)
r2 = jim(axes(ig), fbp_image_r, "FBP image with radon"; clim)
jim(r1,r2)

# For comparison, here is the difference image.
jim(axes(ig), fbp_image_b - fbp_image_r, "Error image")

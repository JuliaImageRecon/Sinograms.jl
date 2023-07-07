#=
# [2D Branchless Distance-driven Projection and Backprojection](@id 08-bdd2d)

This page describes the 2D branchless distance-driven projection
and backprojection method
of
[Basu & De Man, 2006](https://doi.org/10.1117/12.659893)
for fan-beam geometries with a "flat" detector
using the Julia package
[`Sinograms.jl`](https://github.com/JuliaImageRecon/Sinograms.jl).

This page compares the results to the `radon` transform method.
=#

#srcURL

# ### Setup

# Packages needed here.

using ImageGeoms: ImageGeom, fovs, MaskCircle
import ImageGeoms # downsample
using ImagePhantoms: shepp_logan, SheppLoganToft, radon, phantom
using MIRTjim: jim, prompt
import Plots
using Sinograms: project_bdd, backproject_bdd
import Sinograms # downsample, _ar, _dso
using Sinograms: SinoFanFlat, rays, plan_fbp, fbp, Window, Hamming, sino_geom_plot!, angles
using Unitful: mm

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

#=
## Fan-beam sinogram of SheppLoganToft phantom

For illustration,
we start by synthesizing
a fan-beam sinogram
of the `SheppLoganToft` phantom.

For completeness,
we use units (from `Unitful`),
but units are optional.
=#

# Define the sinogram geometry
down = 4 # save time
rg = SinoFanFlat( ; nb = 910, d = 1.0239mm, na = 360, dsd = 949mm, dod = 408mm)
rg = Sinograms.downsample(rg, down)

# Define the image geometry
ig = ImageGeom(MaskCircle(); dims=(512,512), deltas = (1mm,1mm))
ig = ImageGeoms.downsample(ig, down)

# Make a tuple to define imaging geometry for bdd code. todo
geo = (DSD = rg.dsd, DS0 = Sinograms._dso(rg), pSize = ig.deltas[1],
    dSize = rg.d, nPix = ig.dims[1], nDet = rg.nb, angle = Sinograms._ar(rg))

# Examine the geometry to verify the FOV:
pa = jim(axes(ig), ig.mask; prompt=false)
sino_geom_plot!(rg, ig)

#
prompt()

# Ellipse parameters for SheppLoganToft phantom:
μ = 0.01 / mm # typical linear attenuation coefficient
ob = shepp_logan(SheppLoganToft(); fovs = fovs(ig), u = (1, 1, μ));

# Create a phantom image:
testimage = phantom(axes(ig)..., ob);

# Show the true phantom image.
pt = jim(axes(ig), testimage)

# Fan-beam sinograms for phantom:
if !@isdefined(sinogramR)
    @time sinogramR = radon(rays(rg), ob)
    # Ensure that sinogram is not truncated
    @assert all(==(0), sum(abs, sinogramR, dims=2)[[1,end]])
end;

if !@isdefined(sinogramB)
    @time sinogramB = project_bdd(reverse(rot180(testimage'), dims=2), geo)
    sinogramB = sinogramB' # todo
end;

p1 = jim(axes(rg), sinogramB; prompt=false,
 title="bdd_2d sinogram", xlabel="r", ylabel="ϕ")
p2 = jim(axes(rg), sinogramR; prompt=false,
 title="Radon test sinogram", xlabel="r", ylabel="ϕ")
p12 = jim(p1,p2)

# Difference of sinogram using `bdd_2d` versus `radon` method.
pd = jim(axes(rg), sinogramR - sinogramB;
 title="Difference sinogram", xlabel="r", ylabel="ϕ")


#=
## Backprojection with `bdd_2d`
# Note that the back-projected image is not a useful reconstruction.
# See next section for image reconstruction via FBP.
=#
if !@isdefined(imageB)
    @time imageB = backproject_bdd(sinogramB'/1mm, geo) # todo
    imageB = rotr90(imageB) # todo
end
pb = jim(axes(ig), imageB;
 title="bdd_2d backprojection", xlabel="x", ylabel="y")


#=
## Image reconstruction via FBP
Compare the reconstructed image using the `bdd_2d` sinogram
and the `radon` sinogram.

Here we start with a "plan"
that would save work if we were reconstructing many images.
For illustration we include `Hamming` window.
=#

plan = plan_fbp(rg, ig; window = Window(Hamming(), 1.0))
fbp_image_b = fbp(plan, sinogramB)
fbp_image_r = fbp(plan, sinogramR);

# A narrow color window is needed to see the soft tissue structures:
clim = (0.04, 1.1) .* μ
r1 = jim(axes(ig), fbp_image_b, "FBP image with bdd_2d"; clim)
r2 = jim(axes(ig), fbp_image_r, "FBP image with radon"; clim)
r12 = jim(r1,r2; size=(1000,300))

# For comparison, here is the difference image.
pe = jim(axes(ig), fbp_image_b - fbp_image_r, "Difference image")

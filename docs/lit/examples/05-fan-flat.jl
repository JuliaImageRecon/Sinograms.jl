#=
# [Fan-beam tomography: flat detector](@id 05-fan-flat)

This page describes fan-beam tomographic image reconstruction
using the Julia package
[`Sinograms.jl`](https://github.com/JuliaImageRecon/Sinograms.jl).

This page focuses on fan-beam with a "flat" detector,
i.e., one row of a flat-panel detector
as used in many cone-beam CT (CBCT) systems.

This page was generated from a single Julia file:
[05-fan-flat.jl](@__REPO_ROOT_URL__/05-fan-flat.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`05-fan-flat.ipynb`](@__NBVIEWER_ROOT_URL__/05-fan-flat.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`05-fan-flat.ipynb`](@__BINDER_ROOT_URL__/05-fan-flat.ipynb).


# ### Setup

# Packages needed here.

using Plots: plot, gui # these 2 must precede Sinograms for Requires to work!
using Unitful: cm
using Sinograms: SinoFanFlat, rays, plan_fbp, Window, Hamming, fbp, sino_geom_plot!
using ImageGeoms: ImageGeom, fovs, MaskCircle
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using MIRTjim: jim, prompt


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
### Fan-beam sinogram of Shepp-Logan phantom

For illustration,
we start by synthesizing
a fan-beam sinogram
of the Shepp-Logan phantom.

For completeness,
we use units (from Unitful),
but units are optional.

=#

# Use `ImageGeom` to define the image geometry.
ig = ImageGeom(MaskCircle(); dims=(128,126), deltas = (0.2cm,0.2cm) )

# Use `SinoFanFlat` to define the sinogram geometry.
rg = SinoFanFlat( ; nb = 130, d = 0.3cm, na = 100, dsd = 50cm, dod = 14cm)

# Examine the geometry to verify the FOV:
jim(axes(ig), ig.mask; prompt=false)
sino_geom_plot!(rg, ig)

#
prompt()

# Ellipse parameters for Shepp-Logan phantom:
μ = 0.1 / cm # typical linear attenuation coefficient
ob = shepp_logan(SheppLogan(); fovs = fovs(ig), u = (1, 1, μ))

# Arc fan-beam sinogram for Shepp-Logan phantom:
sino = radon(rays(rg), ob)
jim(axes(rg), sino; title="Shepp-Logan sinogram", xlabel="r", ylabel="ϕ")


#=
## Image reconstruction via FBP
Here we start with a "plan",
which would save work if we were reconstructing many images.
For illustration we include `Hamming` window.
=#

plan = plan_fbp(rg, ig; window = Window(Hamming(), 1.0))
fbp_image, sino_filt = fbp(plan, sino)


# A narrow color window is needed to see the soft tissue structures:
clim = (0.9, 1.1) .* μ
jim(axes(ig), fbp_image, "FBP image for flat case"; clim)

# For comparison, here is the ideal phantom image
true_image = phantom(axes(ig)..., ob, 2)
jim(axes(ig)..., true_image, "True phantom image"; clim)

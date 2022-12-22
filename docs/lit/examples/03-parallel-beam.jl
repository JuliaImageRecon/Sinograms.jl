#=
# [Parallel-beam tomography](@id 03-parallel-beam)

This page describes parallel-beam tomographic image reconstruction
using the Julia package
[`Sinograms.jl`](https://github.com/JuliaImageRecon/Sinograms.jl).

This page was generated from a single Julia file:
[03-parallel-beam.jl](@__REPO_ROOT_URL__/03-parallel-beam.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`03-parallel-beam.ipynb`](@__NBVIEWER_ROOT_URL__/03-parallel-beam.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`03-parallel-beam.ipynb`](@__BINDER_ROOT_URL__/03-parallel-beam.ipynb).


# ### Setup

# Packages needed here.

using Sinograms: SinoPar, rays, plan_fbp, fbp
using ImageGeoms: ImageGeom, fovs, MaskCircle
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using Unitful: mm
using MIRTjim: jim, prompt


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
### Parallel-beam sinogram of Shepp-Logan phantom

For illustration,
we start by synthesizing
a parallel-beam sinogram
of the Shepp-Logan phantom.

For completeness,
we use units (from Unitful),
but units are optional.

=#

# Use `ImageGeom` to define the image geometry.
ig = ImageGeom(MaskCircle(); dims=(128,126), deltas = (2mm,2mm) )

# Use `SinoPar` to define the sinogram geometry.
sg = SinoPar( ; nb = 130, d = 2mm, na = 100)

# Ellipse parameters for Shepp-Logan phantom:
μ = 0.01 / mm # typical linear attenuation coefficient
ob = shepp_logan(SheppLogan(); fovs = fovs(ig), u = (1, 1, μ))

# Radon transform of Shepp-Logan phantom:
sino = radon(rays(sg), ob)
jim(axes(sg), sino; title="Shepp-Logan sinogram", xlabel="r", ylabel="ϕ")

## Image reconstruction via FBP
# Here we start with a "plan",
# which would save work if we were reconstructing many images.

plan = plan_fbp(sg, ig)
fbp_image, sino_filt = fbp(plan, sino)


# A narrow color window is needed to see the soft tissue structures:
clim = (0.9, 1.1) .* μ
jim(axes(ig), fbp_image, "FBP image"; clim)

# For comparison, here is the ideal phantom image
true_image = phantom(axes(ig)..., ob, 2)
jim(axes(ig)..., true_image, "True phantom image"; clim)

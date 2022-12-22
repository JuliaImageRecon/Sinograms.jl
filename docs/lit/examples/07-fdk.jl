#=
# [CBCT FDK](@id 07-fdk)

This page describes image reconstruction
for cone-beam computed tomography (CBCT)
with both arc and flat detectors
using the Julia package
[`Sinograms.jl`](https://github.com/JuliaImageRecon/Sinograms.jl).

This page was generated from a single Julia file:
[07-fdk.jl](@__REPO_ROOT_URL__/07-fdk.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`07-fdk.ipynb`](@__NBVIEWER_ROOT_URL__/07-fdk.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`07-fdk.ipynb`](@__BINDER_ROOT_URL__/07-fdk.ipynb).


# ### Setup

# Packages needed here.

using Plots: plot, gui # these 2 must precede Sinograms for Requires to work!
using Unitful: cm
using Sinograms: CtFanArc, CtFanFlat # CtPar
using Sinograms: rays, plan_fbp, Window, Hamming, fdk, ct_geom_plot3
using ImageGeoms: ImageGeom, MaskCircle, fovs
using ImagePhantoms: ellipsoid_parameters, ellipsoid
using ImagePhantoms: radon, phantom
using MIRTjim: jim, prompt


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
### CBCT projections of 3D Shepp-Logan phantom

For illustration,
we start by synthesizing
CBCT projections
of the 3D Shepp-Logan phantom.

For completeness,
we use units (from Unitful),
but units are optional.
=#

# Use `ImageGeom` to define the image geometry.
ig = ImageGeom(MaskCircle(); dims=(64,62,30), deltas = (1,1,2) .* 0.4cm)

# Ellipsoid parameters for 3D Shepp-Logan phantom:
μ = 0.1 / cm # typical linear attenuation coefficient
params = ellipsoid_parameters( ; fovs = fovs(ig), u = (1, 1, μ))
ob = ellipsoid(params)
# A narrow grayscale window is needed to see the soft tissue structures
clim = (0.95, 1.05) .* μ
# Here is the ideal phantom image:
oversample = 3
true_image = phantom(axes(ig)..., ob, oversample)
jim(axes(ig), true_image, "True 3D Shepp-Logan phantom image"; clim)

# Define the system geometry
# (for some explanation use `?CtGeom`):
p = (ns = 130, ds = 0.3cm, nt = 80, dt = 0.4cm, na = 50, dsd = 200cm, dod = 40cm)
cg = CtFanArc( ; p...)

# Examine the geometry to verify the FOV
# (this is more interesting when interacting via other Plot backends):
ct_geom_plot3(cg, ig)

#
prompt()


# CBCT projections
# using `Sinogram.rays` and `ImagePhantoms.radon`:
proj_arc = radon(rays(cg), ob)
jim(cg.s, cg.t, proj_arc ;
    title="Shepp-Logan projections (arc)", xlabel="s", ylabel="t")


#=
There is no "inverse crime" here
because we compute the projection views
using the analytical phantom geometry,
but then reconstruct
on a discrete grid.

## Image reconstruction via FBP / FDK

We start with a "plan",
which would save work if we were reconstructing many images.
For illustration we include `Hamming` window.
=#

plan = plan_fbp(cg, ig; window = Window(Hamming(), 1.0))
fdk_arc = fdk(plan, proj_arc)
jim(axes(ig), fdk_arc, "FDK image (arc)"; clim)

#
err_arc = fdk_arc - true_image
jim(axes(ig), err_arc, "Error image (arc)"; clim = (-1,1) .* (0.05μ))


#=
## Repeat with flat detector geometry
=#

cg = CtFanFlat( ; p...)
proj_flat = radon(rays(cg), ob)
jim(cg.s, cg.t, proj_flat ;
    title="Shepp-Logan projections (flat)", xlabel="s", ylabel="t")

plan = plan_fbp(cg, ig; window = Window(Hamming(), 1.0))
fdk_flat = fdk(plan, proj_flat)
jim(axes(ig), fdk_flat, "FDK image (flat)"; clim)

#
err_flat = fdk_flat - true_image
jim(axes(ig), err_flat, "Error image (flat)"; clim = (-1,1) .* (0.05μ))

#=
As expected for CBCT,
the largest errors are in the end slices.
=#

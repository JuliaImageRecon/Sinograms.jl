#---------------------------------------------------------
# # [Fan-beam tomography: arc detector](@id 04-fan-arc)
#---------------------------------------------------------

#=
This page describes fan-beam tomographic image reconstruction
using the Julia package
[`Sinograms.jl`](https://github.com/JeffFessler/Sinograms.jl).

This page focuses on fan-beam with an "arc" detector,
i.e., 3rd-generation CT systems
where the detector is an arc
having a focal point
at the X-ray source.
This geometry is popular
in part because it facilitates
anti-scatter grids.

This page was generated from a single Julia file:
[04-fan-arc.jl](@__REPO_ROOT_URL__/04-fan-arc.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`04-fan-arc.ipynb`](@__NBVIEWER_ROOT_URL__/04-fan-arc.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`04-fan-arc.ipynb`](@__BINDER_ROOT_URL__/04-fan-arc.ipynb).


# ### Setup

# Packages needed here.

using Plots: plot, gui # these 3 must precede Sinograms for Requires to work!
using Unitful: mm
using UnitfulRecipes
using Sinograms: SinoFanArc, rays, plan_fbp, fbp, sino_geom_plot!
using ImageGeoms: ImageGeom, fovs, MaskCircle
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using MIRTjim: jim, prompt


# The following line is helpful when running this example.jl file as a script;
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
ig = ImageGeom(MaskCircle(); dims=(128,126), deltas = (2mm,2mm) )

# Use `SinoFanArc` to define the sinogram geometry.
sg = SinoFanArc( ; nb = 130, d = 3.2mm, na = 100, dsd = 400mm, dod = 140mm)

# Examine the geometry to verify the FOV:
jim(axes(ig), ig.mask; prompt=false)
sino_geom_plot!(sg ; ig)

#
prompt()

# Ellipse parameters for Shepp-Logan phantom:
μ = 0.01 / mm # typical linear attenuation coefficient
ob = shepp_logan(SheppLogan(); fovs = fovs(ig), u = (1, 1, μ))

# Flat fan-beam sinogram for Shepp-Logan phantom:
sino = radon(ob).(rays(sg)...)
jim(sg.r, sg.ad, sino; title="Shepp-Logan sinogram", xlabel="r", ylabel="ϕ")


#=
## Image reconstruction via FBP
Here we start with a "plan",
which would save work if we were reconstructing many images.
=#

plan = plan_fbp(sg, ig)
fbp_image, sino_filt = fbp(plan, sino)


# A narrow color window is needed to see the soft tissue structures:
clim = (0.9, 1.1) .* μ
jim(axes(ig), fbp_image, "FBP image for arc case"; clim)

# For comparison, here is the ideal phantom image
true_image = phantom(axes(ig)..., ob, 2)
jim(axes(ig)..., true_image, "True phantom image"; clim)

#=
# [Fan-beam tomography: short scan](@id 06-fan-short)

This page describes fan-beam tomographic image reconstruction
using the Julia package
[`Sinograms.jl`](https://github.com/JuliaImageRecon/Sinograms.jl)
for a "short" scan
where the rotation (`orbit`) is 180° plus the fan angle.
This case requires
[Parker weighting](http://doi.org/10.1118/1.595078).

This page focuses on fan-beam with an "arc" detector,
i.e., 3rd-generation CT systems
where the detector is an arc
having a focal point
at the X-ray source.
This geometry is popular
in part because it facilitates
anti-scatter grids.

This page was generated from a single Julia file:
[06-fan-short.jl](@__REPO_ROOT_URL__/06-fan-short.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`06-fan-short.ipynb`](@__NBVIEWER_ROOT_URL__/06-fan-short.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`06-fan-short.ipynb`](@__BINDER_ROOT_URL__/06-fan-short.ipynb).


# ### Setup

# Packages needed here.

using Plots: plot, gui # these 2 must precede Sinograms for Requires to work!
using Unitful: mm
using Sinograms: SinoFanArc, rays, plan_fbp, fbp, sino_geom_plot!
using ImageGeoms: ImageGeom, fovs, MaskCircle
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using MIRTjim: jim, prompt


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
## Fan-beam sinogram of Shepp-Logan phantom

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

#=
Use `SinoFanArc` to define the sinogram geometry,
with the `:short` option for `orbit` to make a short scan.
Note that even though we specify `na = 100`
we end up with `na = 67` views
because of the `:short` option.
=#
rg = SinoFanArc( :short, ;
    nb = 130, d = 3.2mm, na = 100, dsd = 400mm, dod = 140mm,
)

# Examine the geometry to verify the FOV.
# The `na=67` blue dots show the `:short` scan.
jim(axes(ig), ig.mask; prompt=false)
sino_geom_plot!(rg, ig)

#
prompt()

# Ellipse parameters for Shepp-Logan phantom:
μ = 0.01 / mm # typical linear attenuation coefficient
ob = shepp_logan(SheppLogan(); fovs = fovs(ig), u = (1, 1, μ))

# Short arc fan-beam sinogram for Shepp-Logan phantom:
sino = radon(rays(rg), ob)
jim(axes(rg), sino; title="Shepp-Logan 'short' sinogram",
    xlabel="r", ylabel="ϕ")


#=
## Image reconstruction via FBP
Here we start with a "plan",
which would save work if we were reconstructing many images.
=#

plan = plan_fbp(rg, ig)

# Examine Parker weights:
jim(axes(rg), plan.view_weight; title = "Parker weights",
    xlabel="r", ylabel="ϕ")

# Finally perform FBP:
fbp_image = fbp(plan, sino);

# A narrow color window is needed to see the soft tissue structures:
clim = (0.9, 1.1) .* μ
jim(axes(ig), fbp_image, "FBP image for 'short' arc case"; clim)

# For comparison, here is the ideal phantom image:
true_image = phantom(axes(ig)..., ob, 2)
jim(axes(ig)..., true_image, "True phantom image"; clim)

# Here is the difference image.
# Better sampling would reduce the errors.
jim(axes(ig)..., fbp_image - true_image, "Error image")

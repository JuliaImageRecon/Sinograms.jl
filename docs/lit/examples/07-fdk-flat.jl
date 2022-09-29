#---------------------------------------------------------
# # [CBCT FDK: flat detector](@id 07-fdk-flat)
#---------------------------------------------------------

#=
This page describes image reconstruction
for cone-beam computed tomography (CBCT)
with a flat detector
using the Julia package
[`Sinograms.jl`](https://github.com/JeffFessler/Sinograms.jl).

This page was generated from a single Julia file:
[07-fdk-flat.jl](@__REPO_ROOT_URL__/07-fdk-flat.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`07-fdk-flat.ipynb`](@__NBVIEWER_ROOT_URL__/07-fdk-flat.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`07-fdk-flat.ipynb`](@__BINDER_ROOT_URL__/07-fdk-flat.ipynb).


# ### Setup

# Packages needed here.

using Plots: plot, gui # these 3 must precede Sinograms for Requires to work!
using Unitful: cm
using UnitfulRecipes
using Sinograms: CtFanFlat, rays, plan_fbp, Window, Hamming, fdk #, ct_geom_plot!
using Sinograms: CtFanArc
using Sinograms: CtPar
using ImageGeoms: ImageGeom, MaskCircle, fovs
using ImagePhantoms: ellipsoid # todo cut
using ImagePhantoms: SheppLogan, shepp_logan, radon, phantom
using MIRTjim: jim, prompt


# The following line is helpful when running this example.jl file as a script;
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
ig = ImageGeom(MaskCircle(); dims=(128,126,40), deltas = (0.2cm,0.2cm,0.2cm) )
ig = ImageGeom(MaskCircle(); dims=(64,62,30), deltas = (0.4cm,0.4cm,0.8cm) )

# Define the system geometry.
cg = CtFanFlat( ; ns = 130, nt = 80, ds = 0.3cm, na = 50, dsd = 200cm, dod = 40cm)
cg = CtFanArc( ; ns = 130, nt = 80, ds = 0.3cm, na = 50, dsd = 200cm, dod = 40cm)
#cg = CtPar( ; ns = 130, nt = 80, ds = 0.3cm, na = 50)

#src Examine the geometry to verify the FOV:
#src jim(axes(ig), ig.mask; prompt=false)
#src ct_geom_plot!(sg ; ig) # todo

#
#src prompt()

# Ellipse parameters for 3D Shepp-Logan phantom:
μ = 0.1 / cm # typical linear attenuation coefficient
#src ob = shepp_logan(SheppLogan(); fovs = fovs(ig), u = (1, 1, μ))
ob = [
    ellipsoid((3cm,1cm,2cm), (((1,1) .* 0.43 .* fovs(ig)[1])..., 0.30 * fovs(ig)[3])),
    ellipsoid((3cm,2cm,1cm), (4cm,3cm,2cm), (π/6,0))
]
# Here is the ideal phantom image
clim = (0.9, 1.1) .* μ
clim = (0, 2)
true_image = phantom(axes(ig)..., ob, 3)
#jim(axes(ig)[1:2]..., true_image, "True 3D phantom image"; clim)
# todo: jim aspect_ratio

# CBCT projections
i = rays(cg)
proj = [radon(ob)(i...) for i in i]
jim(cg.s, cg.t, proj; title="Shepp-Logan projections", xlabel="s", ylabel="t")


#=
## Image reconstruction via FBP / FDK
We start with a "plan",
which would save work if we were reconstructing many images.
For illustration we include `Hamming` window. 
=#

plan = plan_fbp(cg, ig; window = Window(Hamming(), 1.0))
fdk_image = fdk(plan, proj)

# A narrow color window is needed to see the soft tissue structures: todo
jim(axes(ig), fdk_image, "FDK image"; clim)


#=
throw()
todo clean up and/or cut
include("cylinder_proj.jl")
include("ellipsoid_proj.jl")
using MIRT
using MIRTjim
#using Plots
#using PlotlyJS
include("feldkamp.jl")
down = 4
fover = 2
cg = ct_geom(:arc ; ns = 256, nt = 240, na = 288, ds = 4, dt = 4, down = down, offset_s = 0.25, offset_t = 0, dsd = 949.075, dod = 408.075)
#println("Value of rmax ")
#println(cg.rmax)

fov = 500
nx = 64
ny = 60
nz = 50
dims = (nx, ny, nz)
delta = fov/nx.*(1, -1, 1)
mask2 = ones(nx,ny) .== ones(nx,ny)
mask2[nx,ny] = 0
mask = repeat(mask2, inner = (1,1,1), outer = (1,1,nz))
println("Image geometry")
#=
subtle difference between matlab code: matlab can handle 3D image geometry being
input into ellipse_im even though it's supposed to be 2D, julia cannot
=#
ig2 = MIRT.ImageGeom(dims = dims[1:2], deltas = delta[1:2], offsets = (0,0), mask = mask[:,:,1])
ig3 = MIRT.ImageGeom(dims = dims, deltas = delta, offsets = (0,0,0), mask = mask)

cyl = [20, 10, 10, 150, 150, Inf, 0, 0]
ell = [20 10 10 150 150 980 0 0 0.01; 80 10 10 50 50 30 0 0 0.01; -10 -40 75 40 40 40 0 0 0.01; -10 80 -20 30 30 30 0 0 0.01]
println("Ellipsoid object initialized")

fover = 2
inds = [1,2,4,5,7,8]
tmp = ellipse_im(ig2, reshape(cyl[inds], (1,6)) ; oversample = fover)
#in line 48 of feldkamp_example.m, tmp is all zeros and has two additional parameters
tmp = repeat(tmp, inner = (1,1,1), outer = (1,1,ig3.nz))
xtrue = tmp + ellipsoid_im(ig3, ell ; oversample = fover)
#jim(xtrue)
#prompt()
#plot(PlotlyJS.heatmap(z=xtrue[:,:,24]))
#Plots.heatmap(xtrue[:,:,24])

tmp = cylinder_proj(cg, reshape(cyl, (1,8)), fover)
proj = tmp + ellipsoid_proj(cg, ell, fover)
#Plots.heatmap(proj[:,:,1])

li_hat = proj

xfdk, junk = feldkamp(cg, ig3, li_hat)
errorvector = xfdk - xtrue
println(minimum(errorvector))
println(maximum(errorvector))

cg.plot3(ig3)
=#

#=
Run the code below to look at Feldkamp results for central and offcenter slices
ix = 1:ig3.nx
iy = convert(Int,ceil(ig3.ny/2))
iz = convert(Int,ceil(ig3.nz/2))
Plots.plot(ix, xtrue[ix,iy,iz], xlim=(1,ig3.nx), ylim=(-0.003, 0.023))
Plots.plot!(ix, xfdk[ix, iy, iz])

iz = 1:ig3.nz
ix = convert(Int, 1+floor(ig3.nx/2))
iy = convert(Int, 1+floor(ig3.ny/2))
Plots.plot(iz, xtrue[ix,iy,iz], xlim=(1,ig3.nz), ylim=(-0.003, 0.023))
Plots.plot!(iz, xfdk[ix, iy, iz])
=#

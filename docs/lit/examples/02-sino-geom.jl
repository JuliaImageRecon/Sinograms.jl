#---------------------------------------------------------
# # [2D sinogram geometry](@id 02-sino-geom)
#---------------------------------------------------------

#=
This page describes the 2D sinogram geometries
available in the Julia package
[`Sinograms.jl`](https://github.com/JeffFessler/Sinograms.jl).

This page was generated from a single Julia file:
[02-sino-geom.jl](@__REPO_ROOT_URL__/02-sino-geom.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`02-sino-geom.ipynb`](@__NBVIEWER_ROOT_URL__/02-sino-geom.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`02-sino-geom.ipynb`](@__BINDER_ROOT_URL__/02-sino-geom.ipynb).


# ### Setup

# Packages needed here.

using Unitful: mm, °
using UnitfulRecipes
using Plots # must precede 'using Sinograms' for sino_plot_rays to work
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat, SinoFan
using Sinograms: sino_plot_rays, rays
using MIRTjim: jim, prompt
using InteractiveUtils: versioninfo


# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
### 2D Sinogram geometries

To perform 2D image reconstruction
from a 2D sinogram,
one must first describe
how the rays in the sinogram are sampled.

Mathematically,
the Radon transform ``p(r,\phi)`` of ``f(x,y)``
is defined for all ``r`` and ``ϕ``
but practical systems
record ``p(r,ϕ)``
only for certain finite sets of samples
of those values.


## Parallel-beam geometry

The simplest form of sinogram sampling
is equally spaced samples
in both ``r`` and ``ϕ``.
This sampling is called a parallel-beam geometry.
(Very few practical systems have this sampling pattern,
but it is an easy place to start.)
Both FBP and iterative reconstruction methods
need to know which values
of ``r`` and ``ϕ``
are sampled.

In this package,
`SinoPar`
is the type that describes
a parallel-beam sinogram geometry.

The built-in defaults provide helpful reminders about the usage.
=#

SinoPar()

#=
This package supports units via the
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl)
and using units is recommended
(but not required).

Here is an example of how to specify
a parallel-beam geometry.
Everything is a named keyword argument with sensible default values.

* `orbit` and `orbit_start` must both be unitless (degrees)
  or have the same units (e.g., degrees or radians).
* detector spacing `d` and `strip_width`
  must both be unitless or have the same units (e.g., mm).
* The projection angles ``ϕ`` are equally space and given by
  `orbit_start + (0:(nb-1))/nb * orbit`.
=#

sg = SinoPar( ;
    nb = 64, # number of radial samples ("bins")
    na = 30, # number of angular samples
    d = 2mm, # detector spacing
    offset = 0.25, # quarter detector offset (unitless)
    orbit = 180, # angular range (in degrees)
    orbit_start = 0, # starting angle (in degrees)
    strip_width = 2mm, # detector width
)

#=
The struct `sg` has numerous useful properties;
type `?SinoGeom` to see the full list.

For example,
to access the angular samples in degrees
type `sg.ad`
=#

sg.ad


# The following function visualizes the sampling pattern.

sino_plot_rays(sg; ylims=(0,180), yticks=(0:90:180), widen=true, title="Parallel")

#
prompt()


#=
## Fan-beam CT with an arc detector (3rd generation CT)

For a fan-beam geometry,
the arguments are the same as for `SinoPar`
with the addition of specifying:
* `dsd` distance from source to detector
* `dod` distance from origin (isocenter) to detector
* `dfs` distance from focal point of detector to source
  (0 for a 3rd gen arc detector, and `Inf` for a flat detector)
* `source_offset` for misaligned systems
   where the ray from the source to the detector center
   does not intersect the isocenter.
   Not fully supported; submit an issue if you need this feature.

Here is an example that corresponds to a GE Lightspeed CT system.
These numbers are published in
[IEEE T-MI Oct. 2006, p.1272-1283](http://doi.org/10.1109/TMI.2006.882141).
=#

sg = SinoFanArc( ; nb=888, na=984,
    d=1.0239mm, offset=1.25, dsd=949.075mm, dod=408.075mm)


# Here is a smaller example for plotting the rays.

sg = SinoFanArc( ; nb=64, na=30,
    d=20mm, offset=0.25, dsd=900mm, dod=400mm)
sino_plot_rays(sg; ylims=(-50,400), yticks=(0:180:360), widen=true,
    title="Fan-beam for arc detector")

#
prompt()


#=
## Fan-beam CT with a flat detector

This geometry is the same as the arc detector
except that `dfs=Inf`.
=#

sg = SinoFanFlat( ; nb=64, na=30,
    d=20mm, offset=0.25, dsd=900mm, dod=400mm)


# Here is its sampling plot
sino_plot_rays(sg; ylims=(-50,400), yticks=(0:180:360), widen=true,
    title="Fan-beam for flat detector")

#
prompt()


#=
## Mojette sampling
This is a specialized sampling geometry
that is currently incompletely supported.
=#

sg = SinoMoj( ; nb=60, na=30)

# Here is its sampling plot
sino_plot_rays(sg; ylims=(0,180), yticks=(0:90:180), widen=true,
    title="Mojette sampling")

#=
Here is a diagram that illustrates how the radial spacing
is a function of the projection view angle for the Mojette geometry.
The key aspect here
is that in each image row
the line intersection lengths are identical.
=#

plot(aspect_ratio=1, xlims=(-1,1) .* 3.5, ylims = (-1,1) .* 2.5,
    xlabel="x", ylabel="x", title = "Mojette line integrals")
default(label="")
for y in -2:2
    plot!([-3, 3], [y, y], color=:black)
end
for x in -3:3
    plot!([x, x], [-2, 2], color=:black)
end
plot_ray(r, ϕ) = plot!(
    r*cos(ϕ) .+ [-1, 1] * 4 * sin(ϕ),
    r*sin(ϕ) .+ [1, -1] * 4 * cos(ϕ),
    color = :blue,
)
ia = 4 # pick an angle
r = rays(sg)[1][:,ia] # radial samples
r = r[abs.(r) .< 3]
ϕ = rays(sg)[2][1,ia] # projection angle
plot_ray.(r, ϕ)
gui()

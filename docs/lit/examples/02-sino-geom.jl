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

using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat, rays
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

The built-in defaults are helpful.
=#

SinoPar()


## todo: more details


## Fan-beam CT with an arc detector (3rd generation CT)

SinoFanArc()


## Fan-beam CT with a flat detector

SinoFanFlat()


## Mojette sampling

SinoMoj()

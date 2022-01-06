# fbp2_example.jl

#using Revise
using Sinograms
using MIRT: image_geom, sino_geom
# todo: ImagePhantoms
using MIRT: ellipse_im_params, ellipse_im, ellipse_sino #, :shepplogan
using MIRTjim: jim
using Plots


# Geometry
down = 2 # down-sampling factor for quicker tests
ig = image_geom(nx=512, ny=504, fov=500)
ig = ig.down(down)

sg = sino_geom(:par, nb=1088, na=984, orbit=360, orbit_start=0, d=541/949, offset=0.25,)
sg = sg.down(down)

# Phantom object
ell = ellipse_im_params(ig, :shepplogan)
xtrue = ellipse_im(ig, ell, oversample = 4, hu_scale=1000)
sino = ellipse_sino(sg, ell, oversample = 4) * 1000 # hu_scale

clim = (1 .+ (-1, 1) .* 0.05) .* 1000
p1 = jim(ig.x, ig.y, xtrue, "Phantom"; clim)
p2 = jim(sg.r, sg.ad, sino, "Sinogram", aspect_ratio=:none)

plan = fbp2(sg, ig)
@info "plan complete"
result,sino_filtered = fbp2(plan, sino)
@info "fbp complete"


p3 = jim(ig.x, ig.y, result, "FBP"; clim = clim)
#p3 = jim(ig.x, ig.y, result)
#p4 = jim(sr.r, sg.ad, sino_filtered, "Filtered Sino")
#p4 = plot(sino_filtered, label="")
p4 = jim(ig.x, ig.y, result-xtrue, "error"; clim=(-1,1).*100)

jim(p1,p2,p3,p4)



#=
#original example:

down = 2
ig = image_geom(nx=512, ny=504, fov=500)
ig = ig.down(down)
sg = sino_geom(:par, nb=888, na=984, orbit=180 , d=541/949, offset=0.25)
sg=sg.down(down)

ell = ellipse_im_params(ig, :shepplogan)
xtrue = ellipse_im(ig, ell, oversample = 4, hu_scale=1000)
sino = ellipse_sino(sg, ell, oversample = 4)

clim = (1 .+ (-1, 1) .* 0.05) .* 1000
p1 = jim(ig.x, ig.y, xtrue; clim)
p2 = jim(sg.r, sg.ad, sino, aspect_ratio=:none)

plan=fbp2(sg, ig)
println("plan complete")
result,sino_filtered=fbp2(plan,sino)
println("fbp complete")

p3 = jim(ig.x, ig.y, result; clim)
p4 = jim(real(sino_filtered))
jim(p1,p2,p3,p4)

=#

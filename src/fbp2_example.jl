using Revise
using Sinograms
using MIRT
using MIRTjim: jim

down = 2
ig = image_geom(nx=512, ny=504, fov=500)
ig = ig.down(down)

sg = sino_geom(:par, nb=888, na=984, orbit=180 , d=541/949, offset=0.25)
sg=sg.down(down)

ell = ellipse_im_params(ig, :shepplogan)

xtrue = ellipse_im(ig, ell, oversample = 4, hu_scale=1000)
sino = ellipse_sino(sg, ell, oversample=4)

clim = (1 .+ (-1, 1) .* 0.05) .* 1000
p1 = jim(ig.x, ig.y, xtrue; clim)
p2 = jim(sg.r, sg.ad, sino, aspect_ratio=:none)

plan=fbp2(sg, ig)
println("plan complete")
result,sino_filtered=fbp2(plan,sino)
println("backprojection complete")

p3 = jim(result)
p4 = jim(sino_filtered)
jim(p1,p2,p3,p4)






function fbp2_example_show()
    down = 2
    ig = image_geom(nx=512, ny=504, fov=500)
    ig = ig.down(down)

    sg = sino_geom(:par, nb=888, na=984, orbit=180, strip_width=dr, dr = 541/949, offset=0.25)

    sg=sg.down(down)

    ell = []
    clim = (1 .+ [-1 1] .* 0.05) .* 1000

    xtrue, ell = ellipse_im(ig, ell, oversample = 4, rot=90, hu_scale=1000)
    sino = ellipse_sino(sg, ell, oversample=4)
end


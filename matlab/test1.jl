using MIRT
using MATLAB
using Plots

# setup matlab to use MIRT
if !@isdefined(irtdir)
	ENV["MATLAB_ROOT"] = "/Applications/MATLAB_R2022a.app"

	irtdir = "/Users/soniakim/Documents/matlab"
	tmp = "addpath('$irtdir')"
	eval_string(tmp)
	mat"setup"
end

mat"ig = image_geom('nx', $(ig.nx), 'ny', $(ig.ny), 'fov', $fov, 'dy', 'dx');"

mat"sg = sino_geom('par', 'nb', $(sg.nb), 'na', $(sg.na), 'dr', $(sg.d), 'strip_width', 'd');"

ell = [30. 0 120 100 0 1]
ell = [0. 0 200 100 0 1]
ell = [0. 40 200 200 0 1;
	   60. 0 120 100 0 1]

@mput ell
mat"im_mat = ellipse_im(ig, ell)"

@mget im_mat

fov = 500.
ig = image_geom( ; nx=512, ny=504, fov=fov)
down = 2 # small size for debug
ig = ig.down(down)

sg = sino_geom(:par, nb=888, na=984, offset=0.25, d=541/959)
sg = sg.down(down)

im_julia = ellipse_im(ig, ell)
plot(jim(im_julia, "julia"), jim(im_mat, "matlab"), jim(im_mat - im_julia)) 

#sino_julia = ellipse_sino(sg, ell)

#mat"[sino_mat, pos, ang] = ellipse_sino(sg, ell)"
#mat"sizeof(sino_mat)"
#@mget sino_mat

#plot(jim(sino_julia, "julia"), jim(sino_mat, "matlab"), jim(sino_mat - sino_julia)) 

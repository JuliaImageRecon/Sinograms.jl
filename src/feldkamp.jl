using MAT
using HDF5
include("fdk_filter.jl")
include("cbct_back.jl")

struct ct_geom
	type
	ns
	nt
	na
	down
	nframe
	frame
	orbit_start
	pitch
	source_z0
	units
	user_source_zs
	orbit
	ds
	dt
	offset_s
	offset_t
	dsd
	dso
	dod
	dfs
end

struct im_geom
	nx
	ny
	nz
	dx
	dy
	dz
	offset_x
	offset_y
	offset_z
	offsets
	fov
	zfov
	down
	mask
	is3
	dim
end

function meshgrid(x,y)
    #basically the equivalent of ndgrid in matlab
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X,Y
end

#line 103 of feldkamp.m
function feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl)
    (ns,nt,na) = size(proj)
    ss = ([-(ns-1)/2:(ns-1)/2;].-offset_s).*ds
    tt = ([-(nt-1)/2:(nt-1)/2;].-offset_t).*dt

    ss, tt = meshgrid(ss, tt)
    if isinf(dfs) #flat
        if w1cyl == 1
            ww1 = dso ./ sqrt.(ss.^2 + tt.^2 .+ dsd^2)
        else
            print("Code should not run here")
        end
    elseif dfs == 0 #arc
        if w1cyl == 1 # weighting that is "exact" for cylindrical-like objects
			ww1 = (dso/dsd) * cos.(ss./dsd) ./ sqrt.((tt./dsd).^2 .+1);
		else
			print("Code should not run here")
		end
	else
		print("Error: other configurations not implemeneted")
	end

	ww1 = real.(ww1)
	for ia = 1:na
		#println(ww1)
		#println(size(proj[:,:,ia] .* ww1))
		proj[:,:,ia] = proj[:,:,ia] .* ww1
	end
	return proj
end

function feldkamp_do(proj, cg, ig, ds, dt, offset_s, offset_t, offset_source, dsd, dso, dfs, orbit, orbit_start, mask, nz, dx, dy, dz, offset_xyz, w1cyl, window, ia_skip, extrapolate_t, use_mex, nthread)
	#step 1 compute the weights
	proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl)

	#step 2 filter each projection view
	(ns, nt, na) = size(proj)
	proj = fdk_filter(proj, window, dsd, dfs, ds)

	#step 3 cone beam backprojection
	img = cbct_back(proj, cg, ig)

	return img, proj
end

function feldkamp(cg, ig, proj)
	return feldkamp_do(proj, cg, ig, cg.ds, cg.dt, cg.offset_s, cg.offset_t, 0, cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ig.mask, ig.nz, ig.dx, ig.dy, ig.dz, [ig.offset_x, ig.offset_y, ig.offset_z], 1, "ramp", 1, 0, 0, 0)
end

test = 1
if test == 1

	# dfs = 0
	# ds = 16
	# dsd = 949
	# dso = 541
	# dt = 16
	# offset_s = 0.25
	# offset_t = 0
	# w1cyl = 1
	# proj = ones(6,4,2)
	# proj[3,2,1] = 5
	# proj[4,3,2] = 7
	# window = "ramp"
	#test each function
	#print(fdk_filter(proj, window, dsd, dfs, ds))
	#k = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs, w1cyl)

	#k = fdk_filter(proj, window, dsd, dfs, ds)
	#file = matopen("fdkexampledata.mat")
	#file = matopen("/Users/jasonhu/Documents/julia_files/feldkamp/testfile.mat")
	vars = matread("/Users/jasonhu/Documents/julia_files/feldkamp/fdkexampledata.mat")
	cg1 = vars["cg"]
	ig1 = vars["ig"]
	xtrue = vars["xtrue"]
	a = cg1["data"]
	b = ig1["data"]
	li_hat = vars["li_hat"]
	proj = vars["proj"]

	#construct a ct_geom object with properties of cg for matlab
	cg = ct_geom(a["type"], a["ns"], a["nt"], a["na"], a["down"], a["nframe"], a["frame"], a["orbit_start"], a["pitch"], a["source_z0"], a["units"], a["user_source_zs"], a["orbit"], a["ds"], a["dt"], a["offset_s"], a["offset_t"], a["dsd"], a["dso"], a["dod"], a["dfs"])
	ig = im_geom(b["nx"], b["ny"], b["nz"], b["dx"], b["dy"], b["dz"], b["offset_x"], b["offset_y"], b["offset_z"], b["offsets"], b["fov"], b["zfov"], b["down"], vars["maskor"], b["is3"], b["dim"])

	xfdk, junk = feldkamp(cg, ig, proj)

	error = xfdk - xtrue
	println(minimum(error))
	println(maximum(error))
end

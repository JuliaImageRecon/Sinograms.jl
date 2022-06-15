using MIRT
using MIRTjim: jim
using MATLAB
using Sinograms: bdd_2d
using Plots

# setup matlab to use MIRT
if !@isdefined(irtdir)
	ENV["MATLAB_ROOT"] = "/Applications/matlab"

	irtdir = "/Users/soniakim/Documents/matlab"
	tmp = "addpath('$irtdir')"
	eval_string(tmp)
	mat"setup"
end

deg = 1
geo = (DSD = 100000, DSO = 99700, pSize = 1, dSize = 0.5, nPix = 256, nDet = 1024,
        theta = deg2rad.(0:deg:360-deg), isoX = 0, isoY = 0)

"""
	phantomImg = shepp_logan(pixel, version)
Generates a shepp-logan phantom image for testing of image reconstruction algorithms,
for a pixel size and a version of shepp-logan
"""
phantomImg = shepp_logan(geo.nPix, SheppLoganToft())
sinogramB = projection(phantomImg',geo)
imageB = backprojection(sinogramB,geo)

imFP_julia = sinogramB'

@mput deg, geo, phantomImg, sinogramB

mat"imFP_mat = sinogramB'"

@mget imFP_mat

plot(jim(1:geo.nDet, 1:length(geo.theta), imFP_julia, "julia_FP"; aspect_ratio = :auto, clim=(0,68)), 
	 jim(1:geo.nDet, 1:length(geo.theta), imFP_mat, "matlab_FP"; aspect_ratio = :auto, clim=(0,68)), 
	 jim(imFP_mat - imFP_julia))

imBP_julia = imageB'

@mput deg, geo, phantomImg, imageB

mat"imBP_mat = imageB'"

@mget imBP_mat

plot!(jim(1:geo.nPix, 1:geo.nPix, imBP_julia, "julia_BP"; aspect_ratio = :auto, clim=(30,95)), 
	  jim(1:geo.nPix, 1:geo.nPix, imBP_mat, "matlab_BP"; aspect_ratio = :auto, clim=(30,95)), 
	  jim(imBP_mat - imBP_julia, "diff"))
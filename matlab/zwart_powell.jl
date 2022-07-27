# zwart_powell - compare Matlab and Julia versions

using MATLAB
using Sinograms: zwart_powell
using LazyGrids: ndgrid_array
#using Plots
#using MIRTjim: jim

# setup matlab to use MIRT
if !@isdefined(irtdir)
    ENV["MATLAB_ROOT"] = "/Applications/freeware/matlab"

    irtdir = "/Users/fessler/src/matlab/alg"
    tmp = "addpath('$irtdir')"
    eval_string(tmp)
    mat"setup"
end


r = LinRange(-1, 1, 2001) * 2
phi = LinRange(0, 2π, 361)
jsino = zwart_powell.(r, phi')

rr, pp = ndgrid_array(r, phi) # for matlab...
@mput rr
@mput pp
mat"msino = ir_radon_zwart_powell(pp, rr)"
@mget msino

@assert jsino ≈ msino

#jim(jim(jsino, "julia"), jim(msino, "matlab"), jim(jsino-msino, "diff"))

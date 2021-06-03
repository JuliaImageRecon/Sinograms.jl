
export fbp2


"""
geom = fbp2(sg,ig,[setup options])
image, sino_filt = fbp2(sino, geom, [recon_options])

FBP 2D tomographic image reconstruction for parallel-beam or fan-beam cases,
with either flat or arc detector for fan-beam case.

To use this, you first call it with the sinogram and image geometry.
The routine returns the initialized "geom" structure.  Thereafter, to
to perform FBP reconstruction, you call this routine with that structure
(perhaps numerous times for the same geometry).
See fbp2_example.m for examples.

in (for setup)
- sg::SinoGeom          
- ig::ImageGeom            

options (for setup)
-type = ""                  type of reconstruction 
                            (not all are implemented for fan-beam)
    * "", "std:mex"         call mex file for fast backprojection
    * "std:mat"             slower matlab backprojector
    * "dsc"                 backproject using Gtomo2_dsc with system9
    * "mojette"             use mojette rebinning and Gtomo2_table
    * "table"               use Gtomo2_table with tabulated linear interp
    * "df, pull"            direct Fourier reconstruction with "pull" interp
-window = ""                "" or "hann", or array
				            if array, then use samples [-K/2, K/2)
-nthread = jf("ncore")

out (for setup)
- geom (struct)             initialized structure

in (for recon)
- sino     [nb na *]        sinogram(s) (line integrals)
- geom                      structure from first call

options (for recon)
- window = ""               "" or "hann", or array
                            if array, then use samples [-K/2, K/2)          

out (for recon)
- image         [nx ny *]	reconstructed image(s)
- sino_filt     [nb na *]	filtered sinogram(s)

"""



function fbp2(args...)
    if ~isnumeric(args[1]) #setup
        return fbp2_setup(args[1],args[2],args[3:end])
    else #recon
        out,sino_filt=fbp2_recon(args[1], args[2], args[3:end])
        #=if nargout 
            varargout[1] = sino_filt;
        end =#
    end
end

#
#fbp2_setup()
#
function fbp2_setup(
    sg::SinoGeom, 
    ig::ImageGeom; 
    type::Symbol="",
    extra::Array=[],
    window::String="",
    nthread::Int=jf("ncore")
    )

    
   
# function fbp2_setup_df_pull(sg, ig, arg)

function fbp2_setup_dsc(sg::SinoGeom, ig::ImageGeom, arg::Array)

end

# function fbp2_setup_moj(sg, ig, arg)

function ir_fbp2_par_parker_wt(sg:: SinoGeom)

end

function fbp2_setup_std(sg::SinoGeom, ig::ImageGeom, arg::Array)

end

#function fbp2_setup_tab(nb, na, arg)

function fbp2_recon(sino::Array, geom; window="")

end

function fbp2_recon_std(sino::Array, geom, )

end

#function fbp2_recon_df_pull(sino, st, arg)
    
#function fbp2_recon_moj(sino, moj)

#function fbp2_make_sino_filter_moj(nb, na, dx, orbit, orbit_start, window)

#function fbp2_apply_sino_filter_moj(sino, H)


 
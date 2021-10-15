export fbp2

using MIRT

abstract type FBPplan end

struct NormalPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Union{Symbol,AbstractVector{<:Real}}
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
    #moj::
end

struct DfPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Union{Symbol,AbstractVector{<:Real}}
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
end

struct MojPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Union{Symbol,AbstractVector{<:Real}}
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
    #moj::
end

struct TabPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Union{Symbol,AbstractVector{<:Real}}
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
end

reale = (x) -> (@assert x ≈ real(x); real(x))


"""
    plan = fbp2(sg, ig; how=:normal, window=:none)

FBP 2D tomographic image reconstruction for parallel-beam or fan-beam cases,
with either flat or arc detector for fan-beam case.
    
To use this, you first call it with the sinogram and image geometry.
The routine returns the initialized "plan". Thereafter, to
to perform FBP reconstruction, you call this routine with the plan
(perhaps numerous times for the same geometry).


in 
- `sg::SinoGeom`
- `ig::ImageGeom`

options
- `how::Symbol`             how to reconstruct
    * `:normal`             default 
    * `:mojette`            use mojette rebinning and Gtomo2_table
- `window::Symbol`          e.g. `:hann` (default: `:none`)
-`T::DataType`              type of sino elements (default: `Float32`)

out
- `plan::FBPplan`            initialized plan

"""
function fbp2(
    sg::SinoGeom,
    ig::ImageGeom;
    how::Symbol=:normal,
    window::Symbol=:none,
    T::DataType=Float32,
    #nthread::Int=jf("ncore")
    )
    
    if how === :normal
        plan = fbp2_setup_normal(sg,ig,window,T)
    elseif how === :dsc
        #plan = fbp2_setup_dsc(sg,ig,how,window)
    elseif how === :df
        #plan = fbp2_setup_df(sg,ig,how,window)
    elseif how === :mojette
        #plan = fbp2_setup_moj(sg,ig,how,window)
    elseif how === :table
        #plan = fbp2_setup_tab(sg,ig,how,window)
    else 
        throw("unknown type: $how")
    end

    return plan
end

function fbp2_par_parker_wt(sg::SinoGeom)
    orbit = abs(sg.orbit)
    na = sg.na
    ad = abs(sg.ad-sg.orbit_start)
    nb=sg.nb

    if orbit<180
        @warn("orbit $orbit < 180")
        return sg.ones
    end

    orbit>360 && throw("only 180 <= orbit <= 360 supported for Parker weighting")
    extra = orbit - 180 #extra beyond 180 

    wt = ones(T,na)
    
    ii = ad .< extra
    wt[ii] = abs2.(sin(ad[ii] ./ extra .* pi ./ 2))
    ii = ad .>= 180
    wt[ii] = abs2.(sin((orbit .- ad[ii]) ./ extra .* pi ./ 2))
    wt = wt * orbit / 180 #trick because of the back-projector normalization
    return repeat(wt, nb, 1) #[nb na] sinogram sized 
end

function fbp2_setup_normal(sg::SinoGeom, ig::ImageGeom, window::Symbol, T::DataType)
    
    weight=ones(T,sg.nb,sg.na)
    if sg isa SinoPar
        if abs(sg.orbit) != 180 && abs(sg.orbit) != 360
            weight = fbp2_par_parker_wt(sg)
        end

    elseif sg isa SinoFan
        sg.orbit != 360 && @warn("short-scan fan-beam Parker weighting not done")
        
    elseif sg isa SinoMoj
        if sg.dx == abs(ig.dx)
            plan.moj.G = Gtomo2_table(sg, ig, ["mojette,back1"], nthread=nthread)
        else
            d=sg.dx
            dx=ig.dx
            @warn("mojette sinogram with d=$d vs image with dx=$dx")
        end
    
        plan.moj.H = fbp2_make_sino_filter_moj(sg.nb, sg.na, sg.dx, sg.orbit, sg.orbit_start, window)

    else 
        throw("bad sino type")
    end
    return NormalPlan(sg,ig,window,weight)
    #TODO: plan.moj

end


"""
    image, sino_filt = fbp2(plan, sino)

recon (returns image)

in 
- `plan::FBPplan`
- `sino::AbstractMatrix{<:Number}`
        
out
- `image::AbstractMatrix{<:Number}`       reconstructed image(s)
- `sino_filt::AbstractMatrix{<:Number}`   filtered sinogram(s)

"""
#=
function fbp2(plan::FBPplan, sino::AbstractMatrix{<:Number})

    plan.sg.dim != size(sino) && throw("bad sino size")
    # comments 
    
    if plan isa NormalPlan
        return fbp2_recon_normal(plan, sino)
    elseif plan.how === :df
        #=
        if ~isempty(opt.window), error 'no window for DF', 
        end
	    image = fbp2_recon_df_pull(sino, geom, opt);
	    sino_filt = [];=#
    elseif plan.how === :mojette
        #=
        if ~isempty(opt.window), error 'window only at setup for mojette', end
	    image = fbp2_recon_moj(sino, geom.moj);
	    sino_filt = [];=#
    elseif plan.how === :table
        throw("not done")
    else 
        throw("type bug")
    end

end 
=#

# fbp2_recon_normal
function fbp2(plan::NormalPlan, sino::AbstractMatrix{<:Number})

    plan.sg.dim != size(sino) && throw("bad sino size")

    if plan.sg isa SinoPar
        sino = sino .* plan.parallel_beam_parker_weight
        
	    sino,_,_,_ = fbp2_sino_filter(:flat, sino, ds = plan.sg.dr, window = plan.window)
        
		#image = fbp2_back(plan.sg, plan.ig, sino, do_r_mask=true) 
        image = fbp2_back(plan.sg, plan.ig, sino) 
        
    elseif plan.sg isa SinoFan

        dfs=plan.sg.dfs

        dfs != 0 && ~isinf(dfs) && throw("only arc or flat fan done")
		if isinf(dfs)
			dtype = :flat
		elseif dfs == 0
			dtype = :arc
		else
			throw("bad detector dfs: $dfs")
		end

		sino = fbp2_sino_weight(plan.sg, sino) #todo 
		sino,_,_,_ = fbp2_sino_filter(dtype, sino,
			ds=plan.sg.ds, dsd=plan.sg.dsd,
			window=plan.window)
		image = fbp2_back(plan.sg, plan.ig, sino)
        
    elseif plan.sg isa SinoMoj #TODO (incomplete)
        
        sino = fbp2_apply_sino_filter_moj(sino, plan.moj.H)

        if plan.sg.dx == abs(plan.ig.dx)
            image = plan.moj.G' * sino # backproject
            image = image * (pi / plan.sg.na) # account for "dphi" in integral
        else # revert to conventional pixel driven
            ig = plan.ig
            sg = plan.sg
            arg1 = [uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, sign(ig.dy) * ig.offset_y] # trick: old backproject
            arg2 = [sg.d(1:sg.na), sg.offset, sg.orbit, sg.orbit_start]
            image = jf_mex("back2", arg1[:], arg2[:], int32(arg.nthread), single(sino))
            image = image .* plan.ig.mask
        end
            
    else
        throw("not done")
    end

    return image, sino
end

#fbp2_recon_moj
function fpb2(plan::MojPlan, sino::AbstractMatrix{<:Number})

    moj=plan.moj
    # bin sinogram to mojette sampling, whether it is fan or parallel
    msino = moj.ob_rebin * sino
    msino = fbp2_apply_sino_filter_moj(msino, moj.H) # filter mojette sinogram
    image = moj.G' * msino # backproject
    return image * (pi / size(sino,2)) # account for "dphi" in integral

end

function fbp_make_sino_filter_moj(nb, na, dx, orbit, orbit_start, window)
    ang = deg2rad(orbit_start .+ (0:na-1)./na .* orbit)
    npad = 2^ceil(log2(2*nb-1)) # padded size

    dr = dx * max(abs(cos(ang)), abs(sin(ang)))
    if true
        junk, H = fbp2_sino_filter(:flat, ones(nb,1), ds=1, window=:window, decon1=0)
            
        # trick: ramp filter usually has 1/dr^2 in it, but convolution by sum
        # requires "dr" so we need dr / dr^2 = 1 / dr
        return H * (1 ./ dr) # [npad,na]
    
    else # trick: make sure cutoff frequencies match
        u0 = 1/2/dx # standard cutoff for coarsest sampling
        r = -(npad/2):(npad/2-1) * dr # [nb na]
        h = u0^2 .* (2 .* sinc.(2*u0*r) - sinc.(u0*r).^2)
        h = h .* repeat(dr, [npad 1]) # extra dr for discrete-space convolution
    #	clf, plot(r, h, '.'), keyboard
        H = reale(fft(fftshift(h,1), 1))
    #	H = fbp_apodize(H, nb, window);
        !isempty(window) && throw("window not done yet due to dr")
    end
end


function fbp2_apply_sino_filter_moj(sino, H)

    reale = (x) -> (@assert x ≈ real(x); real(x))

    nb = size(sino,1)
    npad = 2^ceil(log2(2*nb-1)) # padded size
    sinopad = [sino; zeros(npad-nb,size(sino,2))] # padded sinogram
    sino = ifft(reale(fft(sinopad, 1) .* H), 1) 
    sino = sino[1:nb,:]
end
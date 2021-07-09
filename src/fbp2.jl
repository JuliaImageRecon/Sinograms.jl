export fbp2

using MIRT
using Revise

struct FBPplan
    sg::SinoGeom
    ig::ImageGeom
    how::Symbol
    window::Union{Symbol,AbstractVector{<:Real}}
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
end


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

out
- `plan::FBPplan`            initialized plan



"""
function fbp2(
    sg::SinoGeom,
    ig::ImageGeom;
    how::Symbol=:normal,
    window::Symbol=:none
    #nthread::Int=jf("ncore")
    )
    
    if how === :normal
        plan = fbp2_setup_normal(sg,ig,how,window)
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

    wt = ones(na)
    
    ii = ad .< extra
    wt[ii] = sin(ad[ii] ./ extra .* pi ./ 2).^2
    ii = ad .>= 180
    wt[ii] = sin((orbit .- ad[ii]) ./ extra .* pi ./ 2).^2
    wt = wt * orbit / 180 #trick because of the back-projector normalization
    return repeat(wt, nb, 1) #[nb na] sinogram sized 
end

function fbp2_setup_normal(sg::SinoGeom, ig::ImageGeom, how::Symbol, window::Symbol)
    
    weight=ones(sg.na,sg.nb)'
    if sg isa SinoPar
        if abs(sg.orbit) != 180 && abs(sg.orbit) != 360
            weight = fbp2_par_parker_wt(sg)
        end
    
        
    
    elseif sg isa SinoFan
        sg.orbit != 360 && @warn("short-scan fan-beam Parker weighting not done")
        
        #...
    elseif sg isa SinoMoj
        #...
    else 
        throw("bad sino type")
    end
    plan = FBPplan(sg,ig,how,window,weight)
    return plan
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
function fbp2(plan::FBPplan, sino::AbstractMatrix{<:Number})

    plan.sg.dim != size(sino) && throw("bad sino size")
    # comments 
    
    if plan.how === :normal
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

function fbp2_recon_normal(plan::FBPplan, sino::AbstractMatrix{<:Number})
    if plan.sg isa SinoPar
        sino = sino .* plan.parallel_beam_parker_weight
        
	    sino = fbp2_sino_filter(:flat, sino, ds = plan.sg.dr, window = plan.window)

	    
		image = fbp2_back(plan.sg, plan.ig, sino[1]) # single(sino) ?
	    
	    
        
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
		sino = fbp2_sino_weight(plan.sg, sino)
		sino = fbp2_sino_filter(dtype, sino,
			ds=plan.sg.ds, dsd=plan.sg.dsd,
			window=plan.window)
		return fbp2_back_fan(plan.sg, plan.ig, sino)
        
    elseif plan.sg isa SinoMoj
        #=
        sino = fbp2_apply_sino_filter_moj(sino, geom.moj.H);

        if geom.sg.dx == abs(geom.ig.dx)
            image = geom.moj.G' * sino; % backproject
            image = image * (pi / geom.sg.na); % account for "dphi" in integral
        else % revert to conventional pixel driven
            ig = geom.ig;
            sg = geom.sg;
            arg1 = {uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, ...
                sign(ig.dy) * ig.offset_y}; % trick: old backproject
            arg2 = {sg.d(1:sg.na), sg.offset, sg.orbit, sg.orbit_start};
            image = jf_mex('back2', arg1{:}, arg2{:}, ...
                    int32(arg.nthread), single(sino));
            image = image .* geom.ig.mask;
            end
            =#
    else
        throw("not done")
    end

    return image, sino
end


#function fbp2_apply_sino_filter_moj(sino, H)
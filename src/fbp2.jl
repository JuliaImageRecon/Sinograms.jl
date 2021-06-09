export fbp2

struct FBPplan
    sg::SinoGeom
    ig::ImageGeom
    how::Symbol
    window::AbstractVector{<:Real}
    parallel_beam_parker_weight::
    
end


"""
    plan = fbp2(sg, ig; how, window)

FBP 2D tomographic image reconstruction for parallel-beam or fan-beam cases,
with either flat or arc detector for fan-beam case.
    
To use this, you first call it with the sinogram and image geometry.
The routine returns the initialized "plan" structure. Thereafter, to
to perform FBP reconstruction, you call this routine with that structure
(perhaps numerous times for the same geometry).


in 
- `sg::SinoGeom`
- `ig::ImageGeom`

options
- `how::Symbol=:normal`       how to reconstruct
    * :normal               default 
    * :mojette              use mojette rebinning and Gtomo2_table
- `window::Symbol=:none`      e.g. :hann

out
- `plan::FBPplan`            initialized structure

"""
function fbp2(
    sg::SinoGeom,
    ig::ImageGeom;
    how::Symbol=:normal,
    window::Symbol=:none
    #nthread::Int=jf("ncore")
    )
    
    if how === :normal
        if sg isa SinoPar
            if abs(sg.orbit) != 180 && abs(sg.orbit) != 360
                plan.parallel_beam_parker_weight = fbp2_par_parker_wt(sg)
            end
            #...
        elseif sg isa SinoFan
            sg.orbit==360 && throw("short-scan fan-beam Parker weighting not done")
            
            #...
        elseif sg isa SinoMoj
            #...
        else 
            throw("bad sino type")
        end
        
        
    end
    return plan
end


"""
    image, sino_filt=fbp2(plan, sino)

recon (returns image)
in 
- `plan::FBPplan`
- `sino::AbstractMatrix{<:Number}`
        
out
- `image::AbstractMatrix{<:Number}`       reconstructed image(s)
- `sino_filt::AbstractMatrix{<:Number}`   filtered sinogram(s)


"""
function fbp2(plan::FBPplan, sino::AbstractMatrix{<:Number})

    (plan.sg.nb != plan.sg.dim || plan.sg.na != plan.sg.dim) && throw("bad sino size")
    # comments 
    if plan.how === :normal
        if plan.sg isa SinoPar
            #if isvar 
            sino = fbp2_sino_filter(:flat, sino, :ds, plan.sg.dr, :window, window)
            
            if plan.how === :normal
                #...
            end
        elseif plan.sg isa SinoFan
            plan.sg.dfs != 0 && ~isinf(plan.sg.dfs) && throw("only arc or flat fan done")

            
            #...
        end
        
    end

end 

function fbp2_par_parker_wt(sg::SinoGeom)
    orbit = abs(sg.orbit)
    na = sg.na
    ad = abs(sg.ad-sg.orbit_start)
    nb=sg.nb

    if orbit<180
        @warn("orbit $orbit < 180")
        return 1
    end

    orbit>360 && throw("only 180 <= orbit <= 360 supported for Parker weighting")
    
    extra = orbit - 180 #extra beyond 180 

    wt = ones(na)
    #=
    ii = ad < extra
    wt(ii) = sin(ad(ii) / extra * pi / 2).^2
    ii = ad >= 180
    wt(ii) = sin((orbit - ad(ii)) / extra * pi / 2).^2
    =#
    wt = wt * orbit / 180 #trick because of the back-projector normalization
    wt = repeat(wt, nb, 1) #[nb na] sinogram sized

end

 end
#function fbp2_apply_sino_filter_moj(sino, H)
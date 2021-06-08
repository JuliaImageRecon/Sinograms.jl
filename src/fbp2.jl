export fbp2

struct FBPGeom
    sg::SinoGeom
    ig::ImageGeom
    how::Symbol
    window::AbstractVector{<:Real}
end


"""
    geom = fbp2(sg, ig; how, window)
    FBP 2D tomographic image reconstruction for parallel-beam or fan-beam cases,
    with either flat or arc detector for fan-beam case.
    
    To use this, you first call it with the sinogram and image geometry.
    The routine returns the initialized "geom" structure.  Thereafter, to
    to perform FBP reconstruction, you call this routine with that structure
    (perhaps numerous times for the same geometry).

    setup (returns geom)
        in 
        - sg::SinoGeom
        - ig::ImageGeom

        options
        - how::Symbol=:normal       how to reconstruct
                * :normal           default 
                * :mojette          use mojette rebinning and Gtomo2_table
        - window::Symbol=:none      e.g. :hann

        out
        - geom::FBPGeom             initialized structure

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
                geom.parallel_beam_parker_weight = ir_fbp2_par_parker_wt(sg)
            end
            #...
        elseif sg isa SinoFan
            sg.orbit==360 && throw("short-scan fan-beam Parker weighting not done")
            d(x)=convert(Float64,x)
            #...
        elseif sg isa SinoMoj
            #...
        else 
            throw("bad sino type")
        end
        
        
    end
    return geom
end


"""
    image, sino_filt=fbp2(geom, sino; window)

    recon (returns image)
        in 
        - geom::FBPGeom
        - sino::AbstractMatrix{<:Number}
        
        options
        -window::Symbol=:none       e.g. :hann

        out
        - image::AbstractMatrix{<:Number}       reconstructed image(s)
        - sino_filt::AbstractMatrix{<:Number}   filtered sinogram(s)


"""
function fbp2(geom::FBPGeom, sino::AbstractMatrix{<:Number}; window::Symbol=:none)

    (geom.sg.nb != size(sino, 1) || geom.sg.na != size(sino, 2)) && throw("bad sino size")
    # comments 
    if geom.how === :normal
        if geom.sg isa SinoPar
            #if isvar 
            sino = fbp2_sino_filter(:flat, sino, :ds, geom.sg.dr, :window, window)
            
            if geom.how === :normal
                #...
            end
        elseif geom.sg isa SinoFan
            geom.sg.dfs != 0 && ~isinf(geom.sg.dfs) && throw("only arc or flat fan done")

            
            #...
        end
        
    end

end 

function ir_fbp2_par_parker_wt(sg::SinoGeom)
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

    wt = ones(1, na)
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
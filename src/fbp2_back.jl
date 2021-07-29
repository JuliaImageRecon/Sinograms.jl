
export fbp2_back

"""
    img = fbp2_back(sg, ig, sino; ia_skip)

2D backprojection for FBP.

in
- `sg::SinoGeom`                
- `ig::ImageGeom`
- `sino::AbstractMatrix{<:Number}`      sinogram (line integrals), usually ramp filtered

options
- `ia_skip::Int`                        downsample in angle to save time for quick tests (default: 1)

out
- `img::AbstractMatrix{<:Number}`       reconstructed image


"""
function fbp2_back(sg::SinoGeom, ig::ImageGeom, sino::AbstractMatrix{<:Number}; ia_skip::Int=1)

    if sg isa SinoFan
        return fbp2_back_fan(sg, ig, sino, ia_skip=ia_skip)
    elseif sg isa SinoPar
        return fbp2_back(sg, ig, sino, ia_skip)
    else 
        throw("unknown type")
    end
end

function fbp2_back(sg::SinoGeom, ig::ImageGeom, sino::AbstractMatrix{<:Number}, ia_skip::Int)

    # trick: extra zero column saves linear interpolation indexing within loop!

    nb = size(sino,1) # number of radial bins
    nb != sg.nb && throw("nb size") 
    sino=[sino;zeros(eltype(sino),size(sino,2))']

    
    #xc, yc = ndgrid(ig.x, ig.y)
    xc = repeat(ig.x, 1, length(ig.y))
    yc = repeat(ig.y', length(ig.x), 1)
    rr = sqrt.(abs2.(xc) + abs2.(yc)) # [nx ny]

    rmax = ((sg.nb-1)/2-abs(sg.offset)) * sg.d
    #=
    mask = ig.mask
    if do_r_mask
        mask = mask & (rr < rmax);
    end
    xc = xc(mask(:)); % [np] pixels within mask
    yc = yc(mask(:));
    =#

    cang = cos.(sg.ar)
    sang = sin.(sg.ar)

    img = 0
    for ia=1:ia_skip:sg.na
        # ticker(mfilename, ia, sg.na)

        rr = xc .* cang[ia] + yc .* sang[ia] # [np,1]
        rr = rr ./ sg.d .+ sg.w .+ 1 # unitless bin index, +1 because matlab |  NOTE: still +1 in julia?

        # nearest neighbor interpolation:
    #=
    %	ib = round(bb);
    %	if any(ib < 1 | ib > nb), error 'bug', end
    %	% trick: make out-of-sinogram indices point to those extra zeros
    %	ib(ib < 1 | ib > nb) = nb+1;
    %	img = img + sino(ib, ia) ./ L2;
    =#
        # linear interpolation:
        il = floor.(Int64, rr) # left bin
        #=
        if ~do_r_mask
            il = max(il,1);
            il = min(il,nb);
        end
        =#

        #temp
        il = min.(il,nb)
        il = max.(il,1)

    	# (any(il .< 1) || any(il .>= nb)) && throw("bug")

        wr = rr - il # left weight
        wl = 1 .- wr # right weight
        img = img .+ wl .* sino[il, ia] + wr .* sino[il.+1, ia]
    end

    # img = (deg2rad(sg.orbit) / (sg.na/ia_skip)) * embed(img, mask);
    # img = pi / (sg.na/ia_skip) * embed(img, mask) % 2008-10-14

    return pi / (sg.na/ia_skip) * img # NOTE: possible temp


    
    
end
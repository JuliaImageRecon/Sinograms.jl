export fbp2_sino_weight

"""
    sino = fbp2_sino_weight(sg, ig, sino)

Apply sinogram weighting for first step of 2D fan-beam FBP.

in
- `sg::SinoGeom`
- `sino::AbstractMatrix{<:Number}`

out
- `sino::AbstractMatrix{<:Number}`

"""
function fbp2_sino_weight(sg::SinoGeom, sino::AbstractMatrix{<:Number})

    #=
    idim = size(sino);
    sino = reshape(sino, idim(1), idim(2), []);
    =#

    if isinf(sg.dfs)
        sino = fbp2_sino_weight_flat(sino, sg.s, sg.dsd, sg.dso, sg.source_offset)
    elseif sg.dfs == 0
        sino = fbp2_sino_weight_arc(sino, sg.s, sg.dsd, sg.dso, sg.source_offset)
    else
        throw("only flat and arc done")
    end
end

function fbp2_sino_weight_arc(sino, ss, dsd, dso, source_offset)
    na=size(sino,2)
    gam = ss / dsd
    w1 = abs(dso * cos(gam) - source_offset * sin(gam)) / dsd # 1D weighting
    sino = sino .* repeat(w1, [1 na]) 

end


function fbp2_sino_weight_flat(sino, ss, dsd, dso, source_offset)
    na = size(sino,2)
    gam = atan(ss / dsd)
    w1 = abs(dso * cos(gam) - source_offset * sin(gam)) / dsd # 1D weighting
    sino = sino .* repeat(w1, [1 na])
    
end


    


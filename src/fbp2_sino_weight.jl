# fbp2_sino_weight.jl

export fbp2_sino_weight

using MIRT: SinoFanFlat, SinoFanArc

"""
    sino = fbp2_sino_weight(sg, ig, sino)

Apply sinogram weighting for first step of 2D fan-beam FBP.

in
- `sg::SinoGeom`
- `sino::AbstractMatrix{<:Number}`

out
- `sino::AbstractMatrix{<:Number}`

"""
function fbp2_sino_weight(sg::SinoFanArc, sino::AbstractMatrix{<:Number})
    ss = sg.s
    dsd = sg.dsd
    dso = sg.dso
    source_offset = sg.source_offset

    na = size(sino,2)
    gam = ss / dsd

    w1 = @.(abs(dso * cos(gam) - source_offset * sin(gam))) / dsd # 1D weighting
    return sino .* w1 
end


function fbp2_sino_weight(sg::SinoFanFlat, sino::AbstractMatrix{<:Number})
    ss = sg.s
    dsd = sg.dsd
    dso = sg.dso
    source_offset = sg.source_offset

    na = size(sino,2)
    gam = atan.(ss / dsd) # todo: unify with "fun."
    w1 = abs.(dso * cos.(gam) - source_offset * sin.(gam)) / dsd # 1D weighting
    return sino .* w1
end

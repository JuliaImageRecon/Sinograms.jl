# fbp2/sino-weight.jl

export fbp_sino_weight

#using Sinograms: SinoFan, _gamma, _dso


"""
    fbp_sino_weight(rg::SinoFan)
Return 1D sinogram weighting for first step of 2D fan-beam FBP.
"""
function fbp_sino_weight(rg::SinoFan)
    gam = _gamma(rg)
    dso = _dso(rg)
    return @. abs(dso * cos(gam) - rg.source_offset * sin(gam)) / rg.dsd
end

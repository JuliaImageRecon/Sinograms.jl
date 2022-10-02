#=
sys2/footprint.jl
=#

using ImageGeoms: ImageGeom, fovs

export footprint_size

_footprint_di(ig::ImageGeom) = sqrt(sum(abs2, ig.deltas[1:2])) # diagonal size
_footprint_rfov(ig::ImageGeom) = maximum(fovs(ig)[1:2]) / 2
_footprint_ratio(st::SinoGeom, ig::ImageGeom{2}) = _footprint_di(ig) / st.d
_footprint_ratio(st::CtGeom, ig::ImageGeom{3}) = _footprint_di(ig) / st.ds


"""
    footprint_size(st::Union{SinoGeom,CtGeom}, ig::ImageGeom)
Unitless maximum footprint size (in detector pixels).
"""
footprint_size

footprint_size(st::SinoMoj, ig::ImageGeom{2})::Float32 = 3f0 # strip_width=dr
footprint_size(st::SinoPar, ig::ImageGeom{2})::Float32 = _footprint_di(ig) / st.d
footprint_size(st::CtPar, ig::ImageGeom{3})::Float32 = _footprint_di(ig) / st.ds

function footprint_size(st::Union{SinoFanArc,CtFanArc}, ig::ImageGeom)::Float32
    rfov = _footprint_rfov(ig)
    rfov > 0.99 * st.dso && throw("bad dso")
    return _footprint_ratio(st, ig) * st.dsd / (st.dso - rfov)
end

function footprint_size(st::Union{SinoFanFlat,CtFanFlat}, ig::ImageGeom)::Float32
    smax = maximum(abs, st.s)
    rfov = _footprint_rfov(ig)
    rfov > 0.99 * st.dso && throw("bad dso")
    return _footprint_ratio(st, ig) * sqrt(st.dsd^2 + smax^2) / (st.dso - rfov)
end

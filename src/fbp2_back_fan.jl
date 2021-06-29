
export fbp2_back_fan


"""
    img = fbp2_back_fan(sg, ig, sino; ia_skip)

2D backprojection for fan-beam FBP.

in
- `sg::SinoGeom`                
- `ig::ImageGeom`
- `sino::AbstractMatrix{<:Number}`      sinogram (line integrals) 

options
- `ia_skip::Int`                        downsample in angle to save time for quick tests (default: 1)

out
- `img::AbstractMatrix{<:Number}`       reconstructed image

"""
function fbp2_back_fan(sg::SinoGeom, ig::ImageGeom, sino::AbstractMatrix; ia_skip::Int=1)

    sg isa SinoFan || throw("need fan type")

    if sg.dfs == 0
        is_arc=true
    elseif isinf(dsf)
        is_arc=false
    else
        throw("bad dsf")
    end

    return fbp2_back_fan(sino, sg.orbit, sg.orbit_start, 
	sg.dsd, sg.dso, sg.dfs, sg.ds, sg.offset_s, 
	sg.source_offset, 
	ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, 
	is_arc, ig.mask, ia_skip)

end

function fbp2_back_fan(sino::AbstractMatrix{<:Number}, orbit::Union{Symbol,Real}, orbit_start::Real, 
	dsd::RealU, dso::Real, dfs::Real, ds::RealU, offset::Float32, source_offset::Float32, 
	nx::Cint, ny::Cint, dx::Cfloat, dy::Cfloat, offset_x::Cfloat, offset_y::Cfloat,
    is_arc::Bool, mask::Union{Symbol,AbstractArray{Bool}}, ia_skip::Int)

    rmax=[]

    na,nb=size(sino)
    



end


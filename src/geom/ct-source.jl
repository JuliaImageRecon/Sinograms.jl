#=
geom/ct-source.jl
Types for describing X-ray CT source trajectories.
=#

export CtSource
export CtSourceCircle
export CtSourceHelix
#export CtSourceUser


# types

"""
    CtSource
Abstract type for representing X-ray CT source trajectories.
"""
abstract type CtSource end


"""
    CtSourceCircle
Type for representing circular X-ray CT source trajectory.
"""
struct CtSourceCircle <: CtSource end


"""
    CtSourceHelix
Type for representing helical X-ray CT source trajectories
having constant pitch.

# Fields
* `pitch` helix pitch (unitless fraction of zFOV, default 0)
* `source_z0` initial z position of x-ray source (default 0)
  It should have same units as detector pixels etc.
"""
struct CtSourceHelix{Td <: RealU} <: CtSource
    pitch::Float32
    source_z0::Td
end

#=
struct CtSourceUser{Ts <: AbstractVector{<:RealU}} <: CtSource
    source_zs::Ts
end
=#


# constructors

"""
    CtSourceHelix( ; pitch = 0, source_z0 = 0)
Constructor with named keywords
"""
function CtSourceHelix( ; pitch::RealU = 0, source_z0::RealU = 0)
    return CtSourceHelix(Float32(pitch), source_z0)
end


# methods

Base.show(io::IO, ::MIME"text/plain", src::CtSource) = _show(io, MIME("text/plain"), src)


#=
downsample(src::CtSource, down::Int) = src

function downsample(src::CtSourceUser, down::Int)
    source_zs = src.source_zs[1:down:end]
    return CtSourceUser(source_zs)
end
=#

#=
    user-specified angles
    angles = angles[1:down:end]
=#

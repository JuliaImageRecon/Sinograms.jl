#=
geom/util.jl
=#

const Toffset = Float32

"""
    RayGeom{Td,To}
Parent type for `SinoGeom` and `CtGeom`
"""
abstract type RayGeom{Td,To} end


# promote_type version that checks for compatible units
function _promoter(xs::RealU...)
    all(==(oneunit(xs[1])), oneunit.(xs)) || throw("incompatible units")
    T = promote_type(eltype.(xs)...)
    T = promote_type(T, eltype(1f0 * oneunit(T))) # at least Float32
end



"""
    _show(io::IO, ::MIME"text/plain", st::Any)
Informative way to show fields of a struct (composite type).
"""
function _show(io::IO, ::MIME"text/plain", st::Any)
    println(io, "$(typeof(st)) :")
    for f in fieldnames(typeof(st))
        p = getfield(st, f)
        t = (p isa Number) ? _unit_precision(p) : typeof(p)
        println(io, " ", f, "::", t, " ", p)
    end
end


"""
    _shape(x::AbstractArray, dim:Dims)
Reshape `x` to size `dim` with `:` only if needed.
Not type stable!
"""
_shape(x::AbstractArray, dim::Dims) =
    (length(x) == prod(dim)) ? reshape(x, dim) : reshape(x, dim..., :)


# detector centers: d * ((0:nb-1) .- w)
# seems to be needed to help with type inference
function _lin_range(
    d::Td, w::Toffset, n::Int ;
    T::Type{<:Number} = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td <: RealU}
    return d * LinRange(-w, n - 1 - w, n)
end

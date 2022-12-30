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
    T = promote_type(typeof.(xs)...)
    T = promote_type(T, typeof(1f0 * oneunit(T))) # at least Float32
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

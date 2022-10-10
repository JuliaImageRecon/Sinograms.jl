#=
fbp3/ct-geom.jl
CT geometry for 3D tomographic image reconstruction
that describes the sampling characteristics of a CT scan
using the parallel, fan beam arc, or flat fan beam geometry.
2022-05-11, Jason Hu, Jeff Fessler translated from ct_geom.m in MIRT
=#

export dims, downsample, oversample, rays


# Methods common to all types

dims(st::CtGeom) = (st.ns, st.nt, st.na)::NTuple{3,Int}

ct_geom_ws(st::CtGeom) = ((st.ns-1)//2 + st.offset_s)::Toffset
ct_geom_wt(st::CtGeom) = ((st.nt-1)//2 + st.offset_t)::Toffset

#ct_geom_s(st::CtGeom) = st.ds * ((0:st.ns-1) .- st.ws) # can't infer!?
function ct_geom_s(
    st::CtGeom{Td} ;
    T::DataType = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.ds, st.ws, st.ns)
end

function ct_geom_t(
    st::CtGeom{Td} ;
    T::DataType = eltype(oneunit(Td) * one(Toffset)),
)::LinRange{T,Int} where {Td}
    return _lin_range(st.dt, st.wt, st.nt)
end


# down/up sampling

function _downsample(st::CtGeom, down_s::Int, down_t::Int, down_a::Int)
    ns = 4 * max(st.ns ÷ 4down_s, 1)
    nt = 2 * max(st.nt ÷ 2down_t, 1)
    na = max(st.na ÷ down_t, 1)

    out = (ns, nt, st.ds * down_s, st.dt * down_t, st.offset_s, st.offset_t,
         na, st.orbit, st.orbit_start,
    )
    if st isa CtFan
         out = (out..., st.source_offset, st.dsd, st.dod)
    end
    return (out..., st.src)
end


"""
    downsample(ct, down::Int)
    downsample(ct, down::NTuple{3,Int})
Down-sample CT geometry
(for testing with small problems).
"""
downsample(st::CtGeom, down::Int) = downsample(st, (down,down,down))

function downsample(st::G, down::NTuple{3,Real}) where {G <: CtGeom}
    return all(==(1), down) ? st : G(_downsample(st, down...)...)::G
end


function _oversample(st::CtGeom, over::Int)
    return (
        st.ns * over, st.nt * over,
        st.ds / over, st.dt / over,
        st.offset_s * over, st.offset_t * over,
        st.na, st.orbit, st.orbit_start,
    )
end

"""
    oversample(st, over::Int)

Over-sample CT geometry in "radial" dimension.
"""
function oversample(st::G, over::Int) where {G <: CtParallel}
    return (over == 1) ? st : G(_oversample(st, over)..., st.src)::G
end

function oversample(st::G, over::Int) where {G <: CtFan}
    return (over == 1) ? st :
        G(_oversample(st, over)..., st.source_offset, st.dsd, st.dod, st.src)::G
end


# type inference help:
function _rays_type3(Td,To)
    Tϕ = eltype(oneunit(to_radians([oneunit(To)])[1]))
    return Iterators.ProductIterator{Tuple{
        LinRange{eltype(1f0 * oneunit(Td)), Int},
        LinRange{eltype(1f0 * oneunit(Td)), Int},
        LinRange{Tϕ, Int},
        Tϕ,
    }}
end


"""
    i = rays(st::CtGeom)
Return parallel-beam coordinates of all rays for this CT geometry.
Return type of `i` is a `ProductIterator` that makes tuples of the form
`(u, v, ϕ, θ)`.
To make projections call
`p = [fun(c...) for c in i]` where `fun` is `radon(...)`.
"""
function rays(st::CtPar{Td,To})::_rays_type3(Td,To) where {Td,To}
    u = st.s
    v = st.t
    ϕ = st.ar # / oneunit(eltype(st.ar))
    θ = zero(eltype(ϕ))
    i = Iterators.product(u, v, ϕ, θ)
    return i
end


function rays(st::CtFan{Td,To}) where {Td,To}
    st.src isa CtSourceCircle || throw("non-circular not done")
    s = st.s
    t = st.t
    β = st.ar # / oneunit(eltype(st.ar))
    i = Iterators.product(s, t, β)
    if st isa CtFanArc
        fun = stb -> cb_arc_to_par(stb..., st.dso, st.dod)
    else
        fun = stb -> cb_flat_to_par(stb..., st.dso, st.dod)
    end
    return Iterators.map(fun, i)
end



# basic 3d methods


#ct_geom_cone_angle(st::CtParallel) = 0
ct_geom_cone_angle(st::CtFan) = atan((st.nt * st.dt)/2 / st.dsd)

ct_geom_zfov(st::CtParallel) = st.nt * st.dt
ct_geom_zfov(st::CtFan) = st.dso / st.dsd * st.nt * st.dt


_source_dz_per_view(src::CtSource, na, orbit, zfov::Td) where {Td <: RealU} = zero(Td)
function _source_dz_per_view(src::CtSourceHelix, na, orbit, zfov::Td) where {Td <: RealU}
    na_per_360 = na * (360 / orbit)
    return na == 1 ? zero(Td) : src.pitch * zfov / na_per_360
end

ct_geom_source_dz_per_view(st::CtGeom) =
    _source_dz_per_view(st.src, st.na, st.orbit, st.zfov)


ct_geom_source_zs(st::CtGeom{Td,To,<:CtSourceCircle}) where {Td,To} = fill(zero(Td), st.na)
#ct_geom_source_zs(st::CtGeom{Td,To,<:CtSourceUser}) where {Td,To} = st.src.source_zs

function ct_geom_source_zs(st::CtGeom{Td,To,<:CtSourceHelix}) where {Td,To}
    source_dz = ct_geom_source_dz_per_view(st)
    return (0:st.na-1) * source_dz .+ st.src.source_z0
end

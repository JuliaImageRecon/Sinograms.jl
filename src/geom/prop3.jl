#=
geom/prop3.jl
CtGeom properties
=#


const ct_geom_fun_core = (
    (:ws, _ws),
    (:wt, _wt),

    (:s, _s),
    (:t, _t),

    (:xds, _xds),
    (:yds, _yds),

    (:ad, angles),
    (:ar, st -> to_radians(st.ad)),

    (:rfov, _rfov),
#   (:taufun, st -> ((x,y) -> _tau(st,x,y))),

    (:shape, st -> ((x::AbstractArray) -> reshaper(x, dims(st)))),
    (:unitv, st -> ((arg...; kwarg...) -> _unitv(st, arg...; kwarg...))),

    (:zfov, _zfov),
    (:source_dz_per_view, _source_dz_per_view),
    (:source_zs, _source_zs),
)


const ct_geom_fun_fan = (
    ct_geom_fun_core...,
    # fan-specific
    (:dso, st -> st.dsd - st.dod),
    (:dfs, _dfs),
    (:gamma, _gamma),
    (:gamma_s, st -> (ss -> _gamma(st, ss))),
    (:gamma_max, _gamma_max),
    (:gamma_max_abs, _gamma_max_abs),
    (:orbit_short, _orbit_short),
    (:cone_angle, _cone_angle),
)


_props(::CtPar) = ct_geom_fun_core
_props(::CtFan) = ct_geom_fun_fan

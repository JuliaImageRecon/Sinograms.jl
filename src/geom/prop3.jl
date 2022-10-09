#=
geom/prop3.jl
CtGeom properties
=#


const ct_geom_fun_core = (
    (:ws, ct_geom_ws),
    (:wt, ct_geom_wt),

    (:s, ct_geom_s),
    (:t, ct_geom_t),

    (:xds, _geom_xds),
    (:yds, _geom_yds),

    (:ad, angles),
    (:ar, st -> to_radians(st.ad)),

    (:rfov, _geom_rfov),
#   (:taufun, st -> ((x,y) -> ct_geom_tau(st,x,y))),

    (:shape, st -> ((x::AbstractArray) -> reshaper(x, dims(st)))),
    (:unitv, st -> ((arg...; kwarg...) -> _geom_unitv(st, arg...; kwarg...))),

    (:zfov, ct_geom_zfov),
    (:source_dz_per_view, ct_geom_source_dz_per_view),
    (:source_zs, ct_geom_source_zs),
)


const ct_geom_fun_fan = (
    ct_geom_fun_core...,
    # fan-specific
    (:dso, st -> st.dsd - st.dod),
    (:dfs, _geom_dfs),
    (:gamma, _geom_gamma),
    (:gamma_s, st -> (ss -> _geom_gamma_s(st, ss))),
    (:gamma_max, _geom_gamma_max),
    (:gamma_max_abs, _geom_gamma_max_abs),
    (:orbit_short, _orbit_short),
    (:cone_angle, ct_geom_cone_angle),
)


_props(::CtPar) = ct_geom_fun_core
_props(::CtFan) = ct_geom_fun_fan

#=
geom/prop2.jl
SinoGeom properties
=#


const sino_geom_fun_core = (
    (:w, sino_w),

    (:dr, sino_dr), # convenience aka st.d
    (:ds, sino_dr),

    (:r, sino_s),
    (:s, sino_s), # sample locations ('radial')

    (:xds, _geom_xds),
    (:yds, _geom_yds),

    (:ad, angles),
    (:ar, st -> to_radians(st.ad)),

    (:rfov, _geom_rfov),
    (:taufun, st -> ((x,y) -> sino_geom_tau(st,x,y))),

    (:shape, st -> ((x::AbstractArray) -> _shape(x, dims(st)))),
    (:unitv, st -> ((arg...; kwarg...) -> _unitv(st, arg...; kwarg...))),
)


const sino_geom_fun_fan = (
    sino_geom_fun_core...,
    # fan-specific
    (:dso, st -> st.dsd - st.dod),
    (:dfs, _geom_dfs),
    (:gamma, _gamma),
    (:gamma_s, st -> (ss -> _geom_gamma_s(st, ss))),
    (:gamma_max, _gamma_max),
    (:gamma_max_abs, _gamma_max_abs),
    (:orbit_short, _orbit_short),
)


const sino_geom_fun_moj = (
    sino_geom_fun_core...,
    # moj-specific: angular dependent d
    (:d_moj, st -> (ar -> st.d * max(abs(cos(ar)), abs(sin(ar))))),
    (:d_ang, st -> st.d_moj.(st.ar)),
)


_props(::SinoPar) = sino_geom_fun_core
_props(::SinoMoj) = sino_geom_fun_moj
_props(::SinoFan) = sino_geom_fun_fan

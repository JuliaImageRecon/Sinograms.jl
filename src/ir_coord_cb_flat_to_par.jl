function ir_coord_cb_flat_to_par(ss, tt, dso, dod)
    Ds = dso
    Dd = dod
    Dc = Ds + Dd

    uu = Ds * ss ./ sqrt.(ss.^2 .+ Dc^2)
    vv = Ds * tt ./ sqrt.(ss.^2 + tt.^2 .+ Dc^2) .* Dc ./ sqrt.(ss.^2 .+ Dc^2)

    polar = -atan.(tt ./ sqrt.(ss.^2 .+ Dc^2))

    azim = atan.(ss / Dc)
    return uu, vv, azim, polar
end

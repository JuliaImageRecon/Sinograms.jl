#=
fbp/d-angle.jl
=#


#=
# Determine if orbit is roughly a multiple of 360°
function _is360(ar::AbstractVector{<:RealU})
    last = ar[end] - ar[end-1]
    next = ar[end] + last # predict next angle if there were one more
    jump = next - ar[begin]
    jump = mod2pi(jump + π) - π # (-π,π)
    return abs(jump) < minimum(abs, diff(ar))
end
=#


"""
    _angle_weights(ar::AbstractVector{<:RealU})
Angular weighting to # scale projections by dβ (aka dϕ or da)
for Riemann-like integration.

# input
- `ar (na,)` angles in radians.
# output
For now, simply the scalar `(ar[begin+1]-ar[begin]) / n180`
where `n180` is the number of full multiples of 180°.
This simplifies to `π/na`
for typical 180° and 360° scans.
For a fan-beam short scan (180° + fan angle),
it also simplifies to `π/na`,
where the "excess" is handled by Parker weighting.

# Examples

```jldoctest
julia> rg = SinoPar(); Sinograms._angle_weights(Sinograms._ar(rg)), π/rg.na
(0.015707962f0, 0.015707963267948967)
```

```jldoctest
julia> rg = SinoPar(;orbit=360); Sinograms._angle_weights(Sinograms._ar(rg)), π/rg.na
(0.015707962f0, 0.015707963267948967)
```

```jldoctest
julia> rg = SinoFanArc(:short); Sinograms._angle_weights(Sinograms._ar(rg)), π/rg.na
(0.04842342f0, 0.04487989505128276)
```

For now, only equally spaced views are supported,
but this is where unequal spacing would be handled,
and the output would be of size `(na,)`.
"""
function _angle_weights(ar::AbstractVector{<:RealU})
    da = diff(ar) # typically 2π/na for 360° orbit
    all(≈(da[1]), da) || error("Unequal angles not implemented yet")

    # Determine number of full multiples of 180°
    last = ar[end] - ar[end-1]
    next = ar[end] + last # predict next angle if there were one more
    span = next - ar[begin]
    n180 = max(floor(Int, 1.00001abs(span) / π), 1)

    return abs(da[1] / n180)

#=
    abs(ar[end] - ar[begin]) < 2π || error("multi-360 not done")
    if _is360(ar)
        da1 = ar[begin] - ar[end]
        da1 = mod2pi(da1 + π) - π # (-π,π)
        return abs.([da1; da]) / 2
    else
        da1 = da[1]
        return abs.([da1; da])
    end
=#
end

#=
zwart_powell.jl
Translated by Harshit Nanda from ir_radon_zwart_powell.m in 2022.
Matlab author: Seongjin Yoon, Univ. of Michigan, Jan 2015
=#

export zwart_powell

const TOLERANCE = 1e-07


"""
    output = zwart_powell(r, ϕ)
Analytic 2D Radon transform value of Zwart-Powell box spline,
for radial distance `r` (normalized by pixel size) and angle `ϕ` (in radians).
"""
function zwart_powell(r::Real, ϕ::Real)
    s, c = sincos(ϕ)
    return zwart_powell(r, (c, s, c+s, c-s))
end


function zwart_powell(r::Real, z::NTuple{4,T}) where T <: Real
    N = sum(z -> abs(z) > TOLERANCE, z)
    return BoxSp4(r, z..., N) / factorial(N-1)
end


BoxSp0(y::RealU, N::Int) = y > 0 ? y^(N-1) : zero(y)


BoxSp1(y::RealU, z1::T, N::Int) where T <: Real =
    abs(z1) > TOLERANCE ?
        (BoxSp0(y + 0.5 * z1, N) -
         BoxSp0(y - 0.5 * z1, N)) / z1 :
        BoxSp0(y, N)


BoxSp2(y::RealU, z1::T, z2::T, N::Int) where T <: Real =
    abs(z2) > TOLERANCE ?
        (BoxSp1(y + 0.5 * z2, z1, N) -
         BoxSp1(y - 0.5 * z2, z1, N)) / z2 :
        BoxSp1(y, z1, N)


BoxSp3(y::RealU, z1::T, z2::T, z3::T, N::Int) where T <: Real =
    abs(z3) > TOLERANCE ?
        (BoxSp2(y + 0.5 * z3, z1, z2, N) -
         BoxSp2(y - 0.5 * z3, z1, z2, N)) / z3 :
        BoxSp2(y, z1, z2, N)


function BoxSp4(y::RealU, z1::T, z2::T, z3::T, z4::T, N::Int) where T <: Real
    return abs(z4) > TOLERANCE ?
        (BoxSp3(y + 0.5 * z4, z1, z2, z3, N) -
         BoxSp3(y - 0.5 * z4, z1, z2, z3, N)) / z4 :
        BoxSp3(y, z1, z2, z3, N)
end

# fbp/reale.jl

"""
    _reale(x ; rtol = …)
Return the real part of `x`,
but warn if the imaginary part is too large.
By default this uses a 100× larger `rtol` than `isapprox` due to FFT errors.
"""
function _reale(x ; rtol::Real = sqrt(eps(real(eltype(x)))) * 100)
    isapprox(x, real(x); rtol) ||
        @warn("x not real $(sqrt(sum(abs2, x-real(x))/sum(abs2,x))) rtol=$rtol")
    return real(x)
end

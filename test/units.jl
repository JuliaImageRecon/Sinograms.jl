#=
units.jl
Test conversion of arrays from degrees to radians.
This is tricky because of package extensions.
=#

using Sinograms: to_radians, fft_filter, Sinograms
using Test

function _is_sinograms(x)
    return parentmodule(to_radians, (typeof(x),)) == Sinograms
end

function _is_units(x)
    return parentmodule(to_radians, (typeof(x),)) ==
        Base.get_extension(Sinograms, :Sinograms_Units)
end


# first test the fallback version
@testset "units-no" begin
    aa = [180f0] # unitless (like degrees)
    ar = @inferred to_radians(aa)
    @test ar ≈ [π]
    @test ar[1] isa Float32
    @test _is_sinograms(aa)
end


# now test the Unitful version
using Unitful: °, rad

@testset "units-yes" begin
    # https://github.com/PainterQubits/Unitful.jl/issues/375
    ar = rad(180f0°)
    @test ar.val isa Float32

    aa = [180f0] # unitless (like degrees)
    ar = @inferred to_radians(aa)
    @test ar ≈ [π]
    @test ar[1] isa Float32
    @test _is_sinograms(aa)

    aa = [1f0rad] # radians
    ar = @inferred to_radians(aa)
    @test ar ≈ [1rad]
    @test ar[1].val isa Float32
    @test _is_units(aa)

    aa = [180f0°] # degrees
    ar = @inferred to_radians(aa)
    @test ar ≈ [1π * rad]
    @test ar[1].val isa Float32
    @test _is_units(aa)
end


using Unitful: mm
using FFTW: fft, ifft

function _tester(fun, x, args...)
    y = fun(x, args...)
    xu = x * 1mm
    yu = @inferred fun(xu, args...)
    return yu == y * 1mm
end

_test_fft(x, args...) = _tester(fft, x, args...)
_test_ifft(x, args...) = _tester(ifft, x, args...)

@testset "fft-units" begin
    for fun in [fft, ifft]
        @test _tester(fun, rand(3,4))
        @test _tester(fun, rand(3,4), 2)
        @test _tester(fun, ones(Int, 7))

        x = rand(Complex{Float32}, 3, 4)
        @test _tester(fun, x)
        @test _tester(fun, x, 2)
        @test _tester(fun, ones(Complex{Int16}, 7))
    end

    # fft_filter:
    H = fft(rand(8))
    x = rand(Float32, 8, 2)
    @test _tester(fft_filter, x, H)
    x = rand(Complex{Float32}, 8, 2)
    @test _tester(fft_filter, x, H)
end

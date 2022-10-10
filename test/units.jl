#=
units.jl
Test conversion of arrays from degrees to radians.
This is tricky because of Requires.jl
=#

using Sinograms: to_radians
using InteractiveUtils: @which
using Test

function _is_sinograms(x)
    tmp = string(@which to_radians(x))
    check = findall("Sinograms.jl", tmp)
    return length(check) == 1
end

function _is_units(x)
    tmp = string(@which to_radians(x))
    check = findall("units.jl", tmp)
    return length(check) == 1
end


# first test the fallback version
@testset "units-no" begin
    aa = [180f0] # unitless (like degrees)
    ar = @inferred to_radians(aa)
    @test ar ≈ [π]
    @test ar[1] isa Float32
    isinteractive() && @test _is_sinograms(aa)
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
    isinteractive() && @test _is_sinograms(aa)

    aa = [1f0rad] # radians
    ar = @inferred to_radians(aa)
    @test ar ≈ [1rad]
    @test ar[1].val isa Float32
    isinteractive() && @test _is_units(aa)

    aa = [180f0°] # degrees
    ar = @inferred to_radians(aa)
    @test ar ≈ [1π * rad]
    @test ar[1].val isa Float32
    isinteractive() && @test _is_units(aa)
end

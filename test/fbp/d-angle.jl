#=
test/fbp/d-angle.jl
=#

using Sinograms: SinoPar, SinoFanArc
#using Sinograms: _is360
using Sinograms: _ar, _angle_weights
using Test: @test, @test_throws, @testset, @inferred


#=
@testset "_is360" begin
    rg = SinoPar()
    ar = _ar(rg)
    @test !(@inferred _is360(ar))

    rg = SinoPar(; orbit=355)
    ar = _ar(rg)
    @test !(@inferred _is360(ar))

    rg = SinoPar(; orbit=360)
    ar = _ar(rg)
    @test (@inferred _is360(ar))

    rg = SinoPar(; orbit=-720, orbit_start=100)
    ar = _ar(rg)
    @test (@inferred _is360(ar))
end
=#


@testset "d-angle" begin
    rg = SinoPar()
    ar = _ar(rg)
    aw = @inferred _angle_weights(ar)
    @test aw isa Float32
#   @test aw isa Vector
    @test aw ≈ Float32(π/rg.na)
#   @test all(≈(Float32(π/rg.na)), aw) # uniform weighting

    rg = SinoPar(; orbit=360)
    ar = _ar(rg)
    aw = @inferred _angle_weights(ar)
    @test aw isa Float32
    @test aw ≈ Float32(π/rg.na)

    rg = SinoPar(; orbit=-720, orbit_start=100)
    ar = _ar(rg)
    aw = @inferred _angle_weights(ar)
    @test aw isa Float32
    @test aw ≈ Float32(π/rg.na)
#   @test_throws ErrorException (@inferred _angle_weights(ar))

    rg = SinoFanArc(Val(:ge1) ; orbit = :short)
    ar = _ar(rg)
    aw = @inferred _angle_weights(ar)
    @test aw isa Float32
    @test aw ≈ Float32(π/rg.na*rg.orbit/180) # todo: check
end

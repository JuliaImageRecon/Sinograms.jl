#=
fbp-par.jl
Test basic parallel-beam FBP
=#

using Sinograms: fbp
using Test: @test, @testset, @inferred
#using MIRTjim: jim # debug only


@testset "fbp-par" begin
    nr = 64
    dr = 10 / nr
    r = (-(nr-1)/2):((nr-1)/2) * dr
    na = 128
    ϕ = (0:(na-1))/na * π

    proj1 = r -> abs(r) < 1 ? sqrt(1 - r^2) : 0.
    proj2 = (r, ϕ, x, y) -> proj1(r - (x * cos(ϕ) + y * sin(ϕ)))
    sino = proj2.(r, ϕ', 3, 1)
    #jim(r, ϕ, sino)

    tmp = @inferred fbp(sino)
    @test tmp isa Matrix
end

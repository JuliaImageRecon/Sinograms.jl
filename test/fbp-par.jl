#=
test/fbp-par.jl
Test basic parallel-beam FBP
=#

using Sinograms: fbp
using Test: @test, @testset, @inferred
#using MIRTjim: jim; jim(:prompt, true) # debug only

@testset "fbp-par" begin
    nr = 64
    dr = 10 / nr
    r = ((-(nr-1)/2):((nr-1)/2)) * dr
    na = 128
    ϕ = (0:(na-1))/na * π

    proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : 0. # unit disk sinogram
    proj2 = (r, ϕ, x, y) -> proj1(r - (x * cos(ϕ) + y * sin(ϕ))) # shifted
    sino = proj2.(r, ϕ', 3, 1)
#   jim(r, ϕ, sino)

    recon = @inferred fbp(sino; dr)
    @test recon isa Matrix
#   jim(recon; prompt=false, gui=true)

    area_hat = sum(recon) * dr^2
    area_true = π * 1^2
    area_err = abs(area_hat - area_true) / area_true
    @test area_err < 1e-3
end

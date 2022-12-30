#=
test/fbp2/fbp-par.jl
Test basic parallel-beam FBP
=#

using Sinograms: fbp
using Test: @test, @testset, @inferred
using Unitful: mm
#using MIRTjim: jim; jim(:prompt, true) # debug only

@testset "fbp-par" begin
    u = 1mm
    nr = 64
    dr = 10u / nr
    r = ((-(nr-1)/2):((nr-1)/2)) * dr
    na = 128
    ϕ = (0:(na-1))/na * π

    proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : zero(r) # unit disk sinogram
    proj2 = (r, ϕ, x, y, w) -> rad * proj1((r - (x * cos(ϕ) + y * sin(ϕ)))/w) # shifted
    rad = 2u
    sino = proj2.(r, ϕ', 2u, 1u, rad)
    sino = Float32.(sino)
    dr = Float32(dr)
#   jim(r, ϕ, sino)

    recon = @inferred fbp(sino ; dr)
    @test recon isa Matrix{Float32}
#   jim(recon; prompt=false, gui=true)

    area_hat = sum(recon) * dr^2
    area_true = π * rad^2
    area_err = abs(area_hat - area_true) / area_true
#   @show area_true area_hat area_err
    @test area_err < 2e-4
end

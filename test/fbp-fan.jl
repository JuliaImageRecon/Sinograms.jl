#=
test/fbp-fan.jl
Test fan-beam FBP
=#

using Sinograms: plan_fbp, fbp, SinoFanArc, SinoFanFlat, rays
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred
#using MIRTjim: jim; jim(:prompt, true) # debug only

include("helper.jl")

@testset "fbp-fan" begin
    ig = ImageGeom()
    proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : 0. # unit disk sinogram
    proj2 = (r, ϕ, x, y) -> proj1(r - (x * cos(ϕ) + y * sin(ϕ))) # shifted
    for sg in (SinoFanArc(), SinoFanFlat())
        r, ϕ = rays(sg)
        sino = proj2.(r, ϕ, 3, 1)
#       jim(r, ϕ, sino)

        plan = @inferred plan_fbp(sg, ig)
        recon, sino_filt = @NOTinferred fbp(plan, sino) # todo
        @test recon isa Matrix
#       jim(recon) # todo: peak is too low by about "dr"
        sino3 = reshape(sino, size(sino)..., 1)
        recon3 = @NOTinferred fbp(plan, sino3) # todo
        @test recon3[:,:,1] == recon
    end
end

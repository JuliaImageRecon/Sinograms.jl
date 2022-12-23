#=
test/fbp-fan.jl
Test fan-beam FBP
=#

using Sinograms: plan_fbp, fbp, SinoFanArc, SinoFanFlat, rays
using ImageGeoms: ImageGeom, MaskCircle, circle
using Test: @test, @testset, @inferred
using Unitful: mm
#using MIRTjim: jim; jim(:prompt, true)

include("helper.jl")

@testset "fbp-fan" begin
    proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : zero(r) # unit disk sinogram
    proj2 = (r, ϕ, x, y, w) -> proj1((r - (x * cos(ϕ) + y * sin(ϕ)))/w) # shifted
    for geom in (SinoFanArc, SinoFanFlat)
        sg = geom( ; d = 0.1f0mm)
        ig = ImageGeom(MaskCircle() ; deltas = (1,1) .* 0.1f0mm)
        mask = circle(ig ; r = sg.rfov)
        ig = ImageGeom(ig.dims, ig.deltas, ig.offsets, mask)
#       r, ϕ = rays(sg)
        i = rays(sg)
        rad = 2mm
        fun(rϕ) = rad * proj2.(rϕ..., 2mm, 1mm, rad)
#       sino = rad * proj2.(r, ϕ, 2mm, 1mm, rad)
        sino = [fun(i) for i in i]
#       jim(sg.s, sg.ad, sino)

        plan = @inferred plan_fbp(sg, ig)

        recon, sino_filt = @inferred fbp(plan, sino)
        @test recon isa Matrix
#       jim(axes(ig), recon)

        area_hat = sum(recon) * prod(ig.deltas)
        area_true = π * rad^2
        area_err = abs(area_hat - area_true) / area_true
        @test area_err < 4e-3

        sino3 = reshape(sino, size(sino)..., 1)
        recon3 = @inferred fbp(plan, sino3)
        @test recon3[:,:,1] == recon
    end
end

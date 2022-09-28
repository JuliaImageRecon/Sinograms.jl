#=
test/fbp3/fdk.jl
Test basic cone-beam FBP (aka FDK).
=#

using Sinograms: plan_fbp, fdk
using Sinograms: CtFanArc, CtFanFlat
using ImageGeoms: ImageGeom
using ImagePhantoms: radon, ellipsoid, volume
using Test: @test, @testset, @test_broken, @inferred
using Unitful: mm
#using MIRTjim: jim; jim(:prompt, true) # debug only

@testset "fdk" begin
    u = 1mm
    ob = ellipsoid( (2u,3u,4u), (8u,7u,6u) )
    ig = ImageGeom( (30, 28, 26), (1,1,1) .* 1u )

    for geo in (CtFanArc, CtFanFlat)
#   geo = (CtFanArc, CtFanFlat)[1]
        cg = geo( ; ns=64, nt=40, ds=1u, dt=1.2u, na=32)

#=
        proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : zero(r) # unit sphere proj.
        proj2 = (u, v, ϕ, θ, x, y, z, w) -> w * proj1((r - (x * cos(ϕ) + y * sin(ϕ)))/w) # shifted
        rad = 2u
        sino = proj2.(r, ϕ', 2u, 1u, rad)
        jim(r, ϕ, sino)
=#

        i = rays(cg)
        proj = [radon(ob)(i...) for i in i]

        plan = plan_fbp(cg, ig)
#       recon = @inferred fdk(plan, proj)
        recon = fdk(plan, proj) # todo: NOTinferred
        @test recon isa Array{<:Number, 3}
#       jim(recon; prompt=false, gui=true)

        vol_hat = sum(recon) * prod(ig.deltas)
        vol_true = volume(ob)
        vol_err = abs(vol_hat - vol_true) / vol_true
        @test_broken vol_err < 2e-6 # todo: bug
    end
end
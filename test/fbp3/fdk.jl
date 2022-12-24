#=
test/fbp3/fdk.jl
Test basic cone-beam FBP (aka FDK).
=#

using Sinograms: plan_fbp, fdk, rays
using Sinograms: CtFanArc, CtFanFlat
using Sinograms: fdk_weight_cyl
using ImageGeoms: ImageGeom
using ImagePhantoms: radon, ellipsoid, volume
using Test: @test, @testset, @test_broken, @inferred
using Unitful: mm
#using MIRTjim: jim; jim(:prompt, true) # debug only


@testset "fdk" begin
    u = 1mm
    ob = ellipsoid( (2u,3u,4u), (8u,7u,6u) )
    ig = ImageGeom( (30, 28, 26), (1,1,1) .* 1u )

    for geo in (CtFanArc, CtFanFlat), shorts in [(), (:short,)]
        cg = geo(shorts... ; ns=64, nt=40, ds=1u, dt=1.2u, na=32)

        @inferred fdk_weight_cyl(cg)

        proj = radon(rays(cg), [ob])

        plan = @inferred plan_fbp(cg, ig)
        recon = @inferred fdk(plan, proj)
        @test recon isa Array{<:Number, 3}
#       jim(recon; prompt=false, gui=true)

        vol_hat = sum(recon) * prod(ig.deltas)
        vol_true = volume(ob)
        vol_err = abs(vol_hat - vol_true) / vol_true
#       @show vol_err
        @test vol_err < (shorts == () ? 3e-3 : 2e-2)
    end
end

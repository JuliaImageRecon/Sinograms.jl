# test/bdd_2d.jl

using Sinograms: projection, backprojection
using Sinograms: SinoFanFlat, rays, plan_fbp, fbp, Window, Hamming
import Sinograms # downsample, _ar, _dso
using ImageGeoms: ImageGeom, fovs, MaskCircle
import ImageGeoms # downsample
using ImagePhantoms: SheppLogan, shepp_logan, SheppLoganToft, radon, phantom
using Test: @test, @testset, @inferred
using LinearAlgebra: dot
using Unitful: mm

@testset "bdd_2d" begin
    # Define geometry
    down = 4 # faster test
    ig = ImageGeom(MaskCircle(); dims=(512,512), deltas = (1mm,1mm))
    rg = SinoFanFlat( ; nb = 910, d = 1.0239mm, na = 360, dsd = 949mm, dod = 408mm)
    ig = ImageGeoms.downsample(ig, down)
    rg = Sinograms.downsample(rg, down)

    geo = (DSD = rg.dsd, DS0 = Sinograms._dso(rg), pSize = ig.deltas[1],
        dSize = rg.d, nPix = ig.dims[1], nDet = rg.nb, angle = Sinograms._ar(rg))

    # Ellipse parameters for Shepp-Logan phantom
    μ = 0.01 / mm # typical linear attenuation coefficient
    ob = shepp_logan(SheppLoganToft(); fovs = fovs(ig), u = (1, 1, μ))

    testimage = phantom(axes(ig)..., ob)

    # forward projector
    sinogramB = projection(reverse(rot180(testimage'), dims=2),geo)
    sinogram = sinogramB / 1mm # include units

    # check dimensions
    @test sinogramB isa Matrix
    @test size(sinogramB) == (length(geo.angle), geo.nDet)

    # Check radon test
    sinogramR = radon(rays(rg), ob)

    # back projector
    imageB = backprojection(sinogram,geo)
    @test imageB isa Matrix
    @test size(imageB) == (geo.nPix, geo.nPix)

    # adjoint
    x = rand(geo.nPix, geo.nPix) / 1mm
    sinoX = projection(x, geo)

    y = rand(length(geo.angle), geo.nDet) / 1mm
    imY = backprojection(y, geo)

    @show (dot(sinoX,y), dot(x, imY))
    @test isapprox(dot(sinoX,y), dot(x, imY); rtol = 2e-2)
end

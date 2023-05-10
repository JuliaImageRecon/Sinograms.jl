using Sinograms: projection, backprojection, SinoFanFlat, rays, plan_fbp, fbp, Window, Hamming
using ImagePhantoms: SheppLogan, shepp_logan, SheppLoganToft, radon, phantom
using Test: @test, @testset, @inferred
using LinearAlgebra: dot

using ImageGeoms: ImageGeom, fovs, MaskCircle
using MIRTjim: jim, prompt
using Unitful: mm, @u_str

@testset "bdd_2d" begin
    # Define geometry
    deg = 1
    geo = (DSD = 949mm, DS0 = 541mm, pSize = 1mm, dSize = 1.0239mm, 
        nPix = 512, nDet = 888, angle = deg2rad.(0:deg:360-deg))
    
    prompt()
    ig = ImageGeom(MaskCircle(); dims=(512,512), deltas = (1mm,1mm) )
    rg = SinoFanFlat( ; nb = 888, d = 1.0239mm, na = 360, dsd = 949mm, dod = 408mm)
    
    # Ellipse parameters for Shepp-Logan phantom
    μ = 0.01 / mm # typical linear attenuation coefficient
    ob = shepp_logan(SheppLoganToft(); fovs = fovs(ig), u = (1, 1, μ))
    
    testimage = phantom(axes(ig)..., ob)
    jim(reverse(testimage, dims=2)) # show the phantom image

    sinogramB = projection(reverse(rot180(testimage'), dims=2),geo)
    sinogram = sinogramB*u"mm^-1" # include units
    imageB = backprojection(sinogram,geo)

    # Check dimensions
    @test sinogramB isa Matrix
    @test size(sinogramB) == (length(geo.angle), geo.nDet)
    @test imageB isa Matrix
    @test size(imageB) == (geo.nPix, geo.nPix)

    # Check radon test
    sinogramR = radon(rays(rg), ob)
    p1 = jim(1:geo.nDet, 1:length(geo.angle), sinogramB'; title="bdd_2d sinogram", aspect_ratio = :auto, xlabel="r (mm)", ylabel="ϕ")
    p2 = jim(axes(rg), sinogramR; title="Radon Test", xlabel="r", ylabel="ϕ")
    p3 = jim(axes(rg), sinogramR - sinogramB', title="Error image", aspect_ratio = :auto, xlabel="r", ylabel="ϕ")
    jim(p1,p2,p3)
    
    # Image reconstruction via FBP
    plan = plan_fbp(rg, ig; window = Window(Hamming(), 1.0))
    fbp_image_b = fbp(plan, sinogramB')
    fbp_image_r = fbp(plan, sinogramR)
    clim = (0.04, 1.1) .* μ
    r1 = jim(axes(ig), fbp_image_b, "FBP image with bdd_2d"; clim)
    r2 = jim(axes(ig), fbp_image_r, "FBP image with radon"; clim)
    r3 = jim(axes(ig), fbp_image_b - fbp_image_r, "Error image")
    jim(r1,r2,r3)

    # Backprojection of bdd_2d
    jim(1:geo.nPix, 1:geo.nPix, imageB'; title="bdd_2d backprojection", aspect_ratio = :auto, xlabel="mm", ylabel="mm")

    #check adjoint 
    # x = randn(geo.nPix, geo.nPix)
    # sinoX = projection(x, geo)
    
    # y = randn(length(geo.angle), geo.nDet)
    # imY = backprojection(y, geo)
    
    #@show (dot(sinoX,y), dot(x, imY))
    #@test dot(sinoX,y) ≈ dot(x, imY)

end

using Sinograms: projection, backprojection
using ImagePhantoms: shepp_logan, SheppLoganToft, radon #to-do
using Test: @test, @testset, @inferred
using LinearAlgebra: dot

@testset "bdd_2d" begin
    deg = 1
    geo = (DSD = 949, DS0 = 541, pSize = 1*4, dSize = 1.0239*4, 
        nPix = 512÷4, nDet = 888÷4, angle = deg2rad.(0:deg:360-deg))
    phantomImg = shepp_logan(geo.nPix, SheppLoganToft())
    sinogramB = projection(phantomImg',geo)
    imageB = backprojection(sinogramB,geo)

    #check dimensions
    sino = @inferred projection(phantomImg',geo)
    @test sino isa Matrix
    @test size(sinogramB) == (length(geo.angle), geo.nDet)
    image = @inferred backprojection(sinogramB,geo)
    @test image isa Matrix
    @test size(imageB) == (geo.nPix, geo.nPix)
    
    #check adjoint 
    x = randn(geo.nPix, geo.nPix)
    sinoX = projection(x, geo)
    
    y = randn(length(geo.angle), geo.nDet)
    imY = backprojection(y, geo)
    
    @show (dot(sinoX,y), dot(x, imY))
    @test dot(sinoX,y) ≈ dot(x, imY)

    #@test sinogramB == sinogramB_correct
    #@test imageB == imageB_correct

    #=
    using MIRTjim: jim

    p1 = jim(1:geo.nDet, 1:length(geo.theta), sinogramB'; aspect_ratio = :auto, clim=(0,68))
    p2 = jim(1:geo.nPix, 1:geo.nPix, imageB'; aspect_ratio = :auto, clim=(30,95))
    jim(p1,p2)
    =#

end

using Sinograms: bdd_2d
using Test: @test, @testset, @inferred
# using JLD: load

@testset "bdd_2d" begin
    deg = 1
    geo = (DSD = 100000, DSO = 99700, pSize = 1, dSize = 0.5, nPix = 256, nDet = 1024,
          theta = deg2rad.(0:deg:360-deg), isoX = 0, isoY = 0)
    phantomImg = shepp_logan(geo.nPix, SheppLoganToft())
    sinogramB = projection(phantomImg',geo)
    imageB = backprojection(sinogramB,geo)

    # sinogramB_correct = load("sinomat.jld")["sinomat"]
    # imageB_correct = load("imagemat.jld")["imagemat"]

    sino = @inferred projection(phantomImg',geo)
    @test sino isa Matrix
    @test size(sinogramB) == (360, 1024)
    image = @inferred backprojection(sinogramB,geo)
    @test image isa Matrix
    @test size(imageB) == (256, 256)

    # @test sinogramB == sinogramB_correct
    # @test imageB == imageB_correct

    #=
    using MIRTjim: jim

    p1 = jim(1:geo.nDet, 1:length(geo.theta), sinogramB'; aspect_ratio = :auto, clim=(0,68))
    p2 = jim(1:geo.nPix, 1:geo.nPix, imageB'; aspect_ratio = :auto, clim=(30,95))
    jim(p1,p2)
    =#

end
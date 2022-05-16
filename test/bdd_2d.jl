using Sinograms: bdd_2d
using Test: @test, @testset, @inferred

@testset "bdd_2d" begin
    deg = 1
    geo = (DSD = 100000, DSO = 99700, pSize = 1, dSize = 0.5, nPix = 256, nDet = 1024,
          theta = deg2rad.(0:deg:360-deg), isoX = 0, isoY = 0)
    phantomImg = shepp_logan(geo.nPix, SheppLoganToft())
    sinogramB = projection(phantomImg',geo)
    imageB = backprojection(sinogramB,geo)

    sino = @inferred projection(phantomImg',geo)
    @test sino isa Matrix
    @test size(sinogramB) == (360, 1024)
    image = @inferred backprojection(sinogramB,geo)
    @test image isa Matrix
    @test size(imageB) == (256, 256)

    @test sinogramB[1,1] == 0.0
    @test sinogramB[200,400] == 46.580197553727686
    @test sinogramB[192,341] == 61.34138815209777

    @test imageB[1,1] == 31.0730328125125
    @test imageB[100,13] == 45.6778292877918
    @test imageB[254,203] == 39.02675898771454

    #=
    using MIRTjim: jim

    p1 = jim(1:geo.nDet, 1:length(geo.theta), sinogramB'; aspect_ratio = :auto, clim=(0,68))
    p2 = jim(1:geo.nPix, 1:geo.nPix, imageB'; aspect_ratio = :auto, clim=(30,95))
    jim(p1,p2)
    =#

end
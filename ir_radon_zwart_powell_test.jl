using LazyGrids: ndgrid
using LazyGrids: btime, @timeo # not exported; just for timing tests here
using Test: @test, @testset, @test_throws

@testset "Image reconstruction" begin
    #defining variables
    theta = LinRange(0,pi,181)
    r = LinRange(-1,1,101) * 2
    (tt, rr) = ndgrid(theta, r)
    sino = ir_radon_zwart_powell(tt, rr)
    #testing if sino is a matrix
    @test sino isa Matrix
end
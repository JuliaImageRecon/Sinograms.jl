using LazyGrids: ndgrid
using Sinograms: zwart_powell
using Test: @test, @testset, @test_throws, @inferred

@testset "zwart_powell" begin

    theta = LinRange(0, π, 181)
    r = LinRange(-1, 1, 101) * 2

    # Allocating arrays using ndgrid
    (tt, rr) = ndgrid(theta, r)
    sino = @inferred zwart_powell(tt, rr)

    # testing if a matrix is produced
    @test sino isa Matrix

    @test zwart_powell([0], [0])[1] == 0.75
    @test zwart_powell([π/4], [0])[1] ≈ 1/√2
end

#=
    using Plots

    theta = (0:3)/3 * π
    r = LinRange(-1, 1, 101) * 2
    (tt, rr) = ndgrid(theta, r)
    sino = @inferred zwart_powell(tt, rr)

    p = plot()
	for (i,θ) in enumerate(theta)
        plot!(r, sino[i,:], label="θ=$θ")
    end
=#

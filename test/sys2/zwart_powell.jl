# test/sys2/zwart_powell.jl

using Sinograms: zwart_powell
using Test: @test, @testset, @test_throws, @inferred

@testset "zwart_powell" begin
    r = LinRange(-1, 1, 101) * 2
    ϕ = LinRange(0, π, 181)

    myfun(r, ϕ) = zwart_powell.(r, ϕ')
    sino = @inferred myfun(r, ϕ)
    @test sino isa Matrix
    @test (@inferred zwart_powell(0, 0)) == 0.75
    @test (@inferred zwart_powell(0, π/4)) ≈ 1/√2

#=
    using Plots
    using MIRTjim: jim
    p1 = jim(r, ϕ, sino)

    ϕ = (0:3)/3 * π
    sino = zwart_powell.(r, ϕ')

    p2 = plot()
    for (i,θ) in enumerate(ϕ)
        plot!(r, sino[:,i], label="ϕ=$θ")
    end
    plot(p1, p2)
=#
end

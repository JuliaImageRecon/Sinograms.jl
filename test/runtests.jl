# runtests.jl

using Test: @test, @testset, detect_ambiguities
using Sinograms

@testset "Sinograms" begin
    # todo
    # include("todo.jl")

    @test length(detect_ambiguities(Sinograms)) == 0
end

# test/fbp2/filter.jl

using Sinograms: SinoPar, SinoGeom, fbp_filter
using Unitful: mm
using Test: @test, @testset, @test_throws, @inferred

include("../helper.jl")

@testset "fbp_filter" begin
    for d_ in (2, 2f0, 2.0),
            du in (1, 1mm)
        d = d_ * du
        sg = @inferred SinoPar(; d)
        Hk = @inferred fbp_filter(sg)

        fun(sg::SinoGeom{Td}) where Td = Td
        T = eltype(1 / oneunit(fun(sg)))
        @test Hk isa Vector{T}
        @test length(Hk) == 2 * sg.nb
    end
end

# test/fbp2/filter.jl

using Sinograms: SinoPar, SinoGeom, fbp_filter, fbp_sino_filter
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


@testset "fbp_sino_filter" begin
    for d_ in (2, 2f0, 2.0),
            du in (1, 1mm)

        d = d_ * du
        sg = @inferred SinoPar(; d, nb=8, na=7)
        Hk = @inferred fbp_filter(sg)

        sino = @inferred ones(sg)
        sfilt = @inferred fbp_sino_filter(sino, Hk)
        @test sfilt isa Matrix

        s3 = @inferred ones(sg.dim..., 3)
        sfilt = @inferred fbp_sino_filter(s3, Hk)
        @test sfilt isa Array{T,3} where T
    end
end

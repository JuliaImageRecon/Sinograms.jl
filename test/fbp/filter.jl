# test/fbp/filter.jl

using Sinograms: SinoPar, SinoGeom, fbp_filter, fbp_sino_filter, dims
using Unitful: mm
using Test: @test, @testset, @test_throws, @inferred

include("../helper.jl")

@testset "fbp/filter" begin
    for d_ in (2, 2f0, 2.0),
            du in (1, 1mm)
        d = d_ * du
        rg = @inferred SinoPar(; d)
        Hk = @inferred fbp_filter(rg)

        fun(::SinoGeom{Td}) where Td = Td
        T = typeof(one(1f0 * d) / oneunit(fun(rg)))
        @test Hk isa Vector{T}
        @test length(Hk) == 2 * rg.nb
    end
end


@testset "fbp_sino_filter" begin
    for d_ in (2, 2f0, 2.0), du in (1, 1mm)
        d = d_ * du
        rg = @inferred SinoPar(; d, nb=8, na=7)
        Hk = @inferred fbp_filter(rg)

        sino = @inferred ones(rg)
        sfilt = @inferred fbp_sino_filter(sino, Hk)
#       @code_warntype fbp_sino_filter(sino, Hk)
        @test sfilt isa Matrix

        s3 = @inferred ones(dims(rg)..., 3)
        sfilt = @inferred fbp_sino_filter(s3, Hk)
#       @code_warntype fbp_sino_filter(s3, Hk)
        @test sfilt isa Array{T,3} where T
    end
end

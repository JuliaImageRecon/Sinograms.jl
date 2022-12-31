# test/fbp/filter.jl

using Sinograms: SinoPar, SinoGeom, fbp_filter, fbp_sino_filter, dims
using Sinograms: fft_filter
using FFTW: fft
using Unitful: mm
using Test: @test, @testset, @test_throws, @inferred


@testset "fbp/filter" begin
    for d_ in (2, 2f0, 2.0), du in (1, 1mm)
        d = d_ * du
        rg = @inferred SinoPar(; d)
        Hk = @inferred fbp_filter(rg)

        fun(::SinoGeom{Td}) where Td = Td
        T = typeof(one(1f0 * d) / oneunit(fun(rg)))
        @test Hk isa Vector{T}
        @test length(Hk) == 2 * rg.nb
    end
end


@testset "fft_filter" begin
    for d_ in (2, 2f0, 2.0), du in (1, 1mm), factor in (1, 1+0im)
        d = d_ * du
        nb, na = 16, 3
        data = ones(Float32, nb, na) * factor
        Hk = fft(rand(nb))
        @inferred fft_filter(data, Hk)
    end
end


@testset "fbp_sino_filter" begin
    for d_ in (2, 2f0, 2.0), du in (1, 1mm)
        d = d_ * du
        rg = @inferred SinoPar(; d, nb=8, na=7)
        Hk = @inferred fbp_filter(rg)

        sino = @inferred ones(rg)
        sfilt = @inferred fbp_sino_filter(sino, Hk)
        @test sfilt isa AbstractMatrix

        s3 = @inferred ones(dims(rg)..., 3)
        sfilt = @inferred fbp_sino_filter(s3, Hk)
        @test sfilt isa AbstractArray{T,3} where T
    end
end

#=
test/fbp2/fbp-par.jl
Test basic parallel-beam FBP
=#

using Sinograms: plan_fbp, fbp, SinoPar
using ImageGeoms: ImageGeom, MaskCircle, circle
using Test: @test, @testset, @inferred
using Unitful: mm


@testset "fbp-par" begin
    u = 1mm
    nr = 64
    dr = 10u / nr
    r = ((-(nr-1)/2):((nr-1)/2)) * dr
    na = 128
    ϕ = (0:(na-1))/na * π

    proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : zero(r) # unit disk sinogram
    proj2 = (r, ϕ, x, y, w) -> proj1((r - (x * cos(ϕ) + y * sin(ϕ)))/w) # shift
    rad = 2u
    sino = rad * proj2.(r, ϕ', 2u, 1u, rad)
    sino = Float32.(sino)
    dr = Float32(dr)

    recon = @inferred fbp(sino ; dr)
    @test recon isa Matrix{Float32}

    area_hat = sum(recon) * dr^2
    area_true = π * rad^2
    area_err = abs(area_hat - area_true) / area_true
    @test area_err < 2e-4

    rg = SinoPar( ; d = 0.1f0mm)
    i = rays(rg)
    ig = ImageGeom(MaskCircle() ; deltas = (1,1) .* 0.1f0mm)
    mask = circle(ig ; r = _rfov(rg))
    ig = ImageGeom(ig.dims, ig.deltas, ig.offsets, mask)
    fun(rϕ) = rad * proj2.(rϕ..., 2mm, 1mm, rad) # unitless density
    sino = Float32.([fun(i) for i in i]) # mm units
    plan = @inferred plan_fbp(rg, ig)
    recon = @inferred fbp(plan, sino) # unitless

    sino3 = cat(dims=3, sino, 5sino)
    recon = @inferred fbp(plan, sino3)
    area_hat = vec(mapslices(sum, recon, dims=[1,2]) * prod(ig.deltas)) ./ [1,5]
    area_err = abs.(area_hat .- area_true) / area_true
    @test all(<(5e-4), area_err)
end

using Sinograms: Sinograms
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(Sinograms: Sinograms)
end

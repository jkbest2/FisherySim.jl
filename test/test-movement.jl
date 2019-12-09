@testset "Movement models" begin
    @test eqdist_ap isa PopState
    @test eqdist_ap0.P ≈ eqdist_eigvecs
    @test all(eqdist_ap.P .> 0)
    @test all(sum(move.M; dims = 1) .≈ 1)
    @test all(0 .≤ eqdist_ap0.P .≤ 1)
    @test sum(eqdist_ap0) ≈ 1.0
    @test sum(eqdist_ap) ≈ 100.0
end


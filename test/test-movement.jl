@testset "Movement models" begin
    @test all(sum(move.M; dims = 1) .≈ 1)
    @test all(0 .≤ vecstate(eqdist_ap0) .≤ 1)
    @test sum(eqdist_ap0) ≈ 1.0
    @test sum(eqdist_ap) ≈ 100.0
end


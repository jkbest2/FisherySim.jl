@testset "Simulate" begin
    @test p2 ≈ Psums[2]
    @test Psums[1] == sum(P1)
end


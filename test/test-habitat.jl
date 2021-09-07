@testset "Habitat" begin
    @test length(hab1) == 1
    @test length(hab2) == 2

    @test length(hab1[1][300]) == 1
    @test length(hab1[1][25, 25]) == 1
end

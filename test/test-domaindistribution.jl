@testset "Test domain distributions" begin
    @test size(ln_real) == size(Ω)
    @test size(mvln_real) == size(Ω)

    err_msg = "Passed MidpointDisplacement as a type, need an instance e.g. MidpointDisplacement()"
    @test_throws ErrorException(err_msg) DomainDistribution(MidpointDisplacement, Ω)

    @test size(rand(dd1)) == size(Ω)
    @test size(rand(cdd1)) == size(Ω)
    @test size(rand(cdd2)) == size(Ω)
    @test size(rand(bdd1)) == size(Ω)
    @test size(rand(bdd2)) == size(Ω)
end

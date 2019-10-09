@testset "Test domain distributions" begin
    @test size(ln_real) == size(Ω)
    @test size(mvln_real) == size(Ω)
    @test size(matlognorm_real) == (size(Ω)..., T)
end


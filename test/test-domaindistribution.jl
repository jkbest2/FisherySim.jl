@testset "Test domain distributions" begin
    @test size(ln_real) == size(Ω)
    @test size(mvln_real) == size(Ω)
    # @test length(matlognorm_real) == T
    # @test all(size.(matlognorm_real) .== Ref(size(Ω)))
end


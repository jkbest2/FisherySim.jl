@testset "Covariance kernels" begin
    @test expkern isa AbstractCovarianceKernel
    @test matkern isa AbstractCovarianceKernel
    @test expΣ isa PDMat
    @test expΣ isa PDMat

    n = length(Ω)
    @test size(expΣ) == (n, n)
    @test size(matΣ) == (n, n)
    @test all(diag(expΣ.mat) .== σ²)
    @test all(diag(matΣ.mat) .== σ²)
end


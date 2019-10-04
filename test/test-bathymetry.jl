@testset "Bathymetry models" begin
    Σeb = FisherySim.expcov.(Ω.distances, σ², ρ)
    Σes = FisherySim.map_symm(d -> FisherySim.expcov(d, σ², ρ), Ω.distances)

    @test Σeb[7, 7] == σ²
    @test Σes[7, 7] == σ²
    @test Σes[5, 10] == σ² * exp(-Ω.distances[5, 10] / ρ)
    @test FisherySim.expcov(Ω.distances[1, 10], σ², ρ) == σ² * exp(-Ω.distances[1, 10] / ρ)

    @test bathy_exp[5, 16] ≥ 0
end


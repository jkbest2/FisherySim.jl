@testset "Vessels" begin
    ## Need to have range in a container for broadcasting to work correctly
    @test all(rand_t .∈ Ref(1:length(Ω)))
    @test all(pref_t .∈ Ref(1:length(Ω)))
    # @test 
    # @test all(pref_hist.weights[1] .== pref_hist.weights[end]
    @test q_const[(l = 5, t = 0)] == 0.2
    @test q_const[(l = 200, t = 1)] == 0.2
    # @test q_const[1, 1] == 0.2
    # @test q_const[5, 15] == 0.2
    @test q_spat[(l = 500, t = 100)] > 0
    @test q_sptemp[(l = 400, t = 1)] > 0
    # @test q_vary[35] == 1.0

    @test c1 isa Catch
    @test any(P.P .!= 0.005)

    @test vessels(fleet) == fleet.vessels
    @test c1.catch_biomass ≥ 0

    @test sum(getfield.(c2, :catch_biomass)) > 0
    @test sum(P.P) < 1
end


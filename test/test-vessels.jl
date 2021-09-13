@testset "Vessels" begin
    ## Need to have range in a container for broadcasting to work correctly
    @test all(rand_t .∈ Ref(1:length(Ω)))
    @test all(pref_t .∈ Ref(1:length(Ω)))
    # @test 
    # @test all(pref_hist.weights[1] .== pref_hist.weights[end]
    @test q_const[(l = 5, t = 0)] == Catchability(0.2)
    @test q_const[(l = 200, t = 1)] == Catchability(0.2)
    # @test q_const[1, 1] == Catchability(0.2)
    # @test q_const[5, 15] == Catchability(0.2)
    @test q_spat[(l = 500, t = 100)].catchability > 0
    # @test q_sptemp[(l = 400, t = 1)].catchability > 0
    # @test q_vary[35].catchability == 1.0
    @test q_hab[10].catchability == q_hab.base_catchability * hab2[2][10] * 1.2

    @test c1 isa Catch
    @test any(P.P .!= 0.005)

    @test vessels(fleet) == fleet.vessels
    @test c1.catch_biomass ≥ 0

    @test sum(getfield.(c2, :catch_biomass)) > 0
    @test sum(P.P) < 1

    fleet2 = Fleet(fleet.vessels, fleet.total_effort, [1, 2, 3, 4])
    eff2 = FisherySim.order_effort(fleet2)
    @test issorted(eff2)
end


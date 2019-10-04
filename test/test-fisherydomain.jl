@testset "Fishery Domain" begin
    ## Check distance calculation
    loctuples = [(0.0, 0.0) (5.0, 12.0); (3.0, 4.0) (0.0, 3.0)]
    @test all(FisherySim.calculate_distances(loctuples)[:, 1] .== [0.0, 5.0, 13.0, 3.0])

    randsamp = sample(Ω, 400)
    randhist = fit(Histogram, randsamp, closed = :right)
    ## FIXME: Need a test for the unweighted sampling case
    prefsamp = sample(Ω, Weights(1:10_000), 400)
    prefhist = fit(Histogram, prefsamp, closed = :right)
    ## FIXME: Is this the best test to use here?
    @test prefhist.weights[1] < prefhist.weights[end]
end


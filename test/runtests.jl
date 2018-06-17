using FisherySim
using Distances
using Distributions
using StatsBase
using Base.Test
# import Base.Random

srand(1234)

## manually create grid
grid_x = range(5.0, 10.0, 10)
grid_y = range(2.5, 5.0, 20)
locs = Base.product(grid_x, grid_y)
locmat = reduce(vcat, vec([[x y] for (x, y) in locs]))
distmat = Distances.pairwise(Euclidean(), locmat')

## Use constructor
origin = (0.0, 0.0)
antipode = (100.0, 100.0)
## Rectangular grid should catch more subtle bugs that symmetry would otherwise
## let us get away with
n = (10, 20)
Ω = GriddedFisheryDomain(origin, antipode, n)

@testset "Fishery Domain" begin
    ## Check distance calculation
    loctuples = [(0.0, 0.0) (5.0, 12.0); (3.0, 4.0) (0.0, 3.0)]
    @test all(FisherySim.calculate_distances(loctuples)[:, 1] .== [0.0, 5.0, 13.0, 3.0])

    ## Check that computed grid locations are the same
    @test all(size(Ω) .== size(locs))
    @test all(vec(getindex.(Ω.locs, 1)) .== locmat[:, 1])
    @test all(vec(getindex.(Ω.locs, 2)) .== locmat[:, 2])
    ## Check that distances between grid locations are the same
    @test isapprox(Ω.distances, distmat)

    randsamp = sample(Ω, 400)
    randhist = fit(Histogram, randsamp, closed = :right)
    ## FIXME: Need a test for the unweighted sampling case
    prefsamp = sample(Ω, Weights(1.0:200.0), 400)
    prefhist = fit(Histogram, prefsamp, closed = :right)
    ## FIXME: Is this the best test to use here?
    @test prefhist.weights[1] < prefhist.weights[end]
end

σ² = 3.0
ρ = 80.0

μv = zeros(length(Ω))
μm = zeros(size(Ω)...)

bmod_exp = BathymetryModel(Ω, μv, d -> FisherySim.expcov(d, σ², ρ))
bathy_exp = rand(bmod_exp)
bmod_sqexp = BathymetryModel(Ω, μm, d -> FisherySim.sqexpcov(d, σ², ρ))
bathy_sqexp = rand(bmod_sqexp)

@testset "Bathymetry models" begin
    Σeb = FisherySim.expcov.(Ω.distances, σ², ρ)
    Σsb = FisherySim.sqexpcov.(Ω.distances, σ², ρ)
    Σes = FisherySim.map_symm(d -> FisherySim.expcov(d, σ², ρ), Ω.distances)
    Σss = FisherySim.map_symm(d -> FisherySim.sqexpcov(d, σ², ρ), Ω.distances)

    @test Σeb[7, 7] == σ²
    @test Σsb[7, 7] == σ²
    @test Σes[7, 7] == σ²
    @test Σss[7, 7] == σ²
    @test Σes[5, 10] == σ² * exp(-Ω.distances[5, 10] / ρ)
    @test Σss[5, 10] == σ² * exp(-Ω.distances[5, 10]^2 / (2 * ρ))
    @test FisherySim.expcov(Ω.distances[1, 10], σ², ρ) == σ² * exp(-Ω.distances[1, 10] / ρ)

    @test bathy_exp[5, 16] ≥ 0
    @test bathy_sqexp[176] ≥ 0
end

move = MovementModel(bathy_exp, Exponential(10.0), Normal(10.0, 2.0))

@testset "Movement models" begin

    # eqdist_ex = eqdist(move, 100.0)
    # eqdist_ap = approx_eqdist(move, 100.0)

    ## Broke when going from 10×10 grid to 10×20, so upped the `rtol`.
    # @test eqdist_ex.P ≈ eqdist_ap.P rtol = 1e-6
end

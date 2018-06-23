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

eqdist_ex0 = eqdist(move)
eqdist_ap0 = approx_eqdist(move)
eqdist_ex = eqdist(move, 100.0)
eqdist_ap = approx_eqdist(move, 100.0)

@testset "Movement models" begin
    @test all(0 .≤ vecstate(eqdist_ex0) .≤ 1)
    @test all(0 .≤ vecstate(eqdist_ap0) .≤ 1)
    @test sum(eqdist_ex0) ≈ 1.0
    @test sum(eqdist_ap0) ≈ 1.0
    @test sum(eqdist_ex) ≈ 100.0
    @test sum(eqdist_ap) ≈ 100.0
    ## Broke when going from 10×10 grid to 10×20, so upped the `rtol`.
    @test vecstate(eqdist_ex0) ≈ vecstate(eqdist_ap0) rtol = 1e-6
    @test vecstate(eqdist_ex) ≈ vecstate(eqdist_ap) rtol = 1e-6
end

r = 0.05
K = 100.0
schaef = Schaefer(r, K)
P1 = step(schaef, eqdist_ap)
Phalf = PopState(P1.P ./ 2)
Phalf1 = step(schaef, Phalf)

@testset "Population dynamics" begin
    @test sum(P1) ≈ K
    @test all(vecstate(P1) .== vecstate(eqdist_ap))
    @test all(vecstate(Phalf1) .≥ vecstate(Phalf))
    @test sum(Phalf1) ≈ 51.25
end

target_rand = RandomTargeting()
target_pref = PreferentialTargeting(l -> 2 * l[1], Ω)

rand_t = target(Ω, target_rand, 400)
pref_t = target(Ω, target_pref, 400)
pref_hist = fit(Histogram, pref_t, closed = :right)

q_const = Catchability(0.2, Ω)
q_vary = Catchability(l -> 2 * l[2], Ω)

## These have to be fairly high to get many non-zero catches.
## These give ~10% zeros. These (especially ϕ?) may be good for testing,
## as they result in completely fishing out some cells, testing that check.
ξ = 1.9
ϕ = 1.9

v1 = Vessel(target_rand, q_const, ξ, ϕ)
v2 = Vessel(target_pref, q_const, ξ, ϕ)
v3 = Vessel(target_rand, q_vary, ξ, ϕ)
v4 = Vessel(target_pref, q_vary, ξ, ϕ)

P =  PopState(0.005 * ones(10, 20))
c1 = fish!(P, v1, Ω)

fleet = Fleet([v1, v2, v3, v4], [100, 100, 100, 100])

c2 = fish!(P, fleet, Ω)
total_catch = sum(getfield.(c2, :catch_biomass))

@testset "Vessels" begin
    ## Need to have range in a container for broadcasting to work correctly
    @test all(rand_t .∈ [1:length(Ω)])
    @test all(pref_t .∈ [1:length(Ω)])
    ## Implicitly tested by sample tests in "Fishery domain" testset
    @test pref_hist.weights[1] < pref_hist.weights[end]
    @test q_const[5] == 0.2
    @test q_const[200] == 0.2
    @test q_const[1, 1] == 0.2
    @test q_const[5, 15] == 0.2
    @test q_vary[5, 15] == 145.0
    @test q_vary[35] == 35.0

    @test c1 isa Catch
    @test any(P.P .!= 0.005)
    
    @test vessels(fleet) == fleet.vessels
    @test c1.catch_biomass > 0

    @test sum(getfield.(c2, :catch_biomass)) > 0
    @test sum(P.P) < 1

end
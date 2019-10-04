using FisherySim
using Distances
using Distributions
using StatsBase
using Test
using Random

Random.seed!(1234)

## Use constructor
origin = (0.0, 0.0)
antipode = (100.0, 100.0)
n = (100, 100)
Ω = GriddedFisheryDomain(origin, antipode, n)

@testset "Fishery Domain" begin
    ## Check distance calculation
    loctuples = [(0.0, 0.0) (5.0, 12.0); (3.0, 4.0) (0.0, 3.0)]
    @test all(FisherySim.calculate_distances(loctuples)[:, 1] .== [0.0, 5.0, 13.0, 3.0])

    randsamp = sample(Ω, 400)
    randhist = fit(Histogram, randsamp, closed = :right)
    ## FIXME: Need a test for the unweighted sampling case
    prefsamp = sample(Ω, Weights(randn(10_000)), 400)
    prefhist = fit(Histogram, prefsamp, closed = :right)
    ## FIXME: Is this the best test to use here?
    @test prefhist.weights[1] < prefhist.weights[end]
end

σ² = 3.0
ρ = 40.0

μv = zeros(length(Ω))
μm = zeros(size(Ω)...)

bmod_exp = BathymetryModel(Ω, μv, d -> FisherySim.expcov(d, σ², ρ))
bathy_exp = rand(bmod_exp)

@testset "Bathymetry models" begin
    Σeb = FisherySim.expcov.(Ω.distances, σ², ρ)
    Σes = FisherySim.map_symm(d -> FisherySim.expcov(d, σ², ρ), Ω.distances)

    @test Σeb[7, 7] == σ²
    @test Σes[7, 7] == σ²
    @test Σes[5, 10] == σ² * exp(-Ω.distances[5, 10] / ρ)
    @test FisherySim.expcov(Ω.distances[1, 10], σ², ρ) == σ² * exp(-Ω.distances[1, 10] / ρ)

    @test bathy_exp[5, 16] ≥ 0
end

move = MovementModel(bathy_exp, Exponential(10.0), Normal(10.0, 2.0))

eqdist_ap0 = approx_eqdist(move)
eqdist_ap = approx_eqdist(move, 100.0)

@testset "Movement models" begin
    @test all(sum(move.M; dims = 1) .≈ 1)
    @test all(0 .≤ vecstate(eqdist_ap0) .≤ 1)
    @test sum(eqdist_ap0) ≈ 1.0
    @test sum(eqdist_ap) ≈ 100.0
end

r = 0.05
K = 100.0

schaef = Schaefer(r, K)
P1 = schaef(eqdist_ap)
Phalf = PopState(P1.P ./ 2)
Phalf1 = schaef(Phalf)

pt = PellaTomlinson(r, K, 3.39)
P1_pt = pt(eqdist_ap)
P1half_pt = pt(Phalf1)

unispt = StochasticProduction(PellaTomlinson(r, K, 3.39),
                              LogNormal(-0.2^2 / 2, 0.2))
multispt = StochasticProduction(
    PellaTomlinson(r, K, 3.39),
    FisherySim.mean1MvLogNormal(bmod_exp.D.Σ))
P1_uspt = unispt(P1)
P1_mspt = multispt(P1)

@testset "Population dynamics" begin
    @test sum(P1) ≈ K
    @test all(vecstate(P1) .== vecstate(eqdist_ap))
    @test all(vecstate(Phalf1) .≥ vecstate(Phalf))
    @test sum(Phalf1) ≈ 51.25

    @test P1 isa PopState
    @test P1_pt isa PopState
    @test P1_uspt isa PopState
    @test P1_mspt isa PopState
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

P = PopState(ones(size(Ω)...) / length(Ω))
c1 = fish!(P, v1, Ω)

fleet = Fleet([v1, v2, v3, v4], [100, 100, 100, 100])

c2 = fish!(P, fleet, Ω)
total_catch = sum(getfield.(c2, :catch_biomass))

@testset "Vessels" begin
    ## Need to have range in a container for broadcasting to work correctly
    @test all(rand_t .∈ Ref(1:length(Ω)))
    @test all(pref_t .∈ Ref(1:length(Ω)))
    ## Implicitly tested by sample tests in "Fishery domain" testset
    @test pref_hist.weights[1] < pref_hist.weights[end]
    @test q_const[5] == 0.2
    @test q_const[200] == 0.2
    @test q_const[1, 1] == 0.2
    @test q_const[5, 15] == 0.2
    @test q_vary[5, 15] == 29.0
    @test q_vary[35] == 1.0

    @test c1 isa Catch
    @test any(P.P .!= 0.005)

    @test vessels(fleet) == fleet.vessels
    @test c1.catch_biomass ≥ 0

    @test sum(getfield.(c2, :catch_biomass)) > 0
    @test sum(P.P) < 1
end

p1 = PopState(1e-2 * ones(100, 100))
@show sum(P1)
Pvec, Cvec = simulate(P1, fleet, move, schaef, Ω, 10)

Psums = sum.(Pvec)

function getcatch(CV::Vector{C}) where {C<:Catch}
    timevec = getfield.(CV, :time)
    catchvec = getfield.(CV, :catch_biomass)
    catchtot = zeros(length(unique(timevec)))
    for (t, c) in zip(timevec, catchvec)
        catchtot[t] += c
    end
    catchtot
end

Csums = getcatch(Cvec)

p2 =  begin
          rempop = Psums[1] - Csums[1]
          rempop + rempop * r * (1 - rempop / K)
      end

@testset "Simulate" begin
    @test p2 ≈ Psums[2]
    @test Psums[1] == sum(P1)
end

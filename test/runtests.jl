# Load required package
using FisherySim
using Distances
using Distributions
using StatsBase
using Test
using Random
using PDMats
using LinearAlgebra

## Make reproducible
Random.seed!(1234)

## -----------------------------------------------------------------------------
## Construct a GriddedFisheryDomain
origin = (0.0, 0.0)
antipode = (100.0, 100.0)
n = (100, 100)
Ω = GriddedFisheryDomain(origin, antipode, n)
T = 10 # Number of years

include("test-fisherydomain.jl")

## -----------------------------------------------------------------------------
## Construct covariance kernels and covariance matrices over fishery domains
σ² = 3.0
ρ = 40.0

expkern = ExpCov(σ², ρ)
matkern = Matérn32Cov(σ², ρ)
ar1kern = AR1(1.0, 0.5)

expΣ = cov(expkern, Ω)
matΣ = cov(matkern, Ω)
ar1Σ = cov(ar1kern, T)

include("test-covkernels.jl")

## -----------------------------------------------------------------------------
## Define distributions over the domain
lognorm = DomainDistribution(LogNormal(-0.2^2 / 2, 0.2), Ω)
mvlognorm = DomainDistribution(MvLogNormal, ones(length(Ω)), matΣ, Ω)
matlognorm = DomainDistribution(MatrixLogNormal, ones(length(Ω), T), matΣ, ar1Σ, Ω)

ln_real = rand(lognorm)
mvln_real = rand(mvlognorm)
matlognorm_real = rand(matlognorm)

include("test-domaindistribution.jl")

## -----------------------------------------------------------------------------
## Generate bathymetry
μv = zeros(length(Ω))
μm = zeros(size(Ω)...)

bmod_mat = BathymetryModel(Ω, μv, matkern)
bathy_mat = rand(bmod_mat)

include("test-bathymetry.jl")


## -----------------------------------------------------------------------------
## Construct a movement model
move = MovementModel(bathy_mat, Exponential(10.0), Normal(10.0, 2.0))

eqdist_ap0 = approx_eqdist(move)
eqdist_ap = approx_eqdist(move, 100.0)

include("test-movement.jl")

## -----------------------------------------------------------------------------
## Define Schaefer and PellaTomlinson dynamics models
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
    mvlognorm)
P1_uspt = unispt(P1)
P1_mspt = multispt(P1)

include("test-pop_dynamics.jl")

## -----------------------------------------------------------------------------
## Define vessels and fleets
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

include("test-vessels.jl")

## -----------------------------------------------------------------------------
## Fish a population using above definitions
p1 = PopState(1e-2 * ones(100, 100))
Pvec, Cvec = simulate(P1, fleet, move, schaef, Ω, T)

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

p2 = begin
         rempop = Psums[1] - Csums[1]
         rempop + rempop * r * (1 - rempop / K)
     end

include("test-simulate.jl")


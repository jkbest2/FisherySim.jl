using FisherySim
using Distances
using Distributions
using Base.Test

srand(1234)

grid = range(5.0, 10.0, 10)
locs = Base.product(grid, grid)
locmat = reduce(vcat, vec([[x y] for (x, y) in locs]))
distmat = Distances.pairwise(Euclidean(), locmat')

σ² = 3.0
ρ = 80.0

@test FisherySim.expcov(distmat[1, 1], σ², ρ) == σ²
@test FisherySim.expcov(distmat[1, 10], σ², ρ) == σ² * exp(-distmat[1, 10] / ρ)

μ = ones(size(grid, 1)^2)
Σ = FisherySim.expcov.(distmat, σ², ρ)
# write your own tests here
@test all(diag(Σ) .== σ²)
@test Σ[1, 10] == σ² * exp(-distmat[1, 10] / ρ)

bathy_model = BathymetryModel(MvNormal(μ, Σ), locmat, [10.0, 10.0])
bathy = rand(bathy_model)

@test bathy.bathy[1] == 4.190834994994973

move = MovementModel(bathy, Exponential(10.0), Normal(10.0, 2.0))

eqdist_ex = eqdist(move)
eqdist_ap = PopState(approx_eqdist(move))

@test eqdist_ex.P ≈ eqdist_ap.P

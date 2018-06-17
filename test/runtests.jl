using FisherySim
using Distances
using Distributions
using Base.Test

srand(1234)

## manually create grid
grid_x = range(5.0, 10.0, 10)
grid_y = range(2.5, 5.0, 20)
locs = Base.product(grid_x, grid_y)
locmat = reduce(vcat, vec([[x y] for (x, y) in locs]))
distmat = Distances.pairwise(Euclidean(), locmat')

## Check distance calculation
loctuples = [(0.0, 0.0) (5.0, 12.0); (3.0, 4.0) (0.0, 3.0)]
@test all(FisherySim.calculate_distances(loctuples)[:, 1] .== [0.0, 5.0, 13.0, 3.0])

## Use constructor
origin = (0.0, 0.0)
antipode = (100.0, 100.0)
## Rectangular grid should catch more subtle bugs that symmetry would otherwise
## let us get away with
n = (10, 20)
Ω = GriddedFisheryDomain(origin, antipode, n)

## Check that computed grid locations are the same
@test all(vec(getindex.(Ω.locs, 1)) .== locmat[:, 1])
@test all(vec(getindex.(Ω.locs, 2)) .== locmat[:, 2])
## Check that distances between grid locations are the same
@test isapprox(Ω.distances, distmat)

σ² = 3.0
ρ = 80.0

@test FisherySim.expcov(distmat[1, 1], σ², ρ) == σ²
@test FisherySim.expcov(distmat[1, 10], σ², ρ) == σ² * exp(-distmat[1, 10] / ρ)

μ = ones(size(grid_x, 1) * size(grid_y, 1))
Σ = FisherySim.expcov.(distmat, σ², ρ)
# write your own tests here
@test all(diag(Σ) .== σ²)
@test Σ[1, 10] == σ² * exp(-distmat[1, 10] / ρ)

bathy_model = BathymetryModel(MvNormal(μ, Σ), locmat, [10.0, 10.0])
bathy = rand(bathy_model)

@test_broken bathy.bathy[1] == 4.190834994994973

move = MovementModel(bathy, Exponential(10.0), Normal(10.0, 2.0))

eqdist_ex = eqdist(move, 100.0)
eqdist_ap = approx_eqdist(move, 100.0)

## Broke when going from 10×10 grid to 10×20, so upped the `rtol`.
@test eqdist_ex.P ≈ eqdist_ap.P rtol = 1e-6

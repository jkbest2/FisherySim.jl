"""
    MovementModel
        M::Matrix{Float64}

Model for movement preferences.
"""
struct MovementModel
    M::Matrix{Float64}
end

function MovementModel(B::Bathymetry, distance::Distributions.UnivariateDistribution)
    distmat = pairwise(Euclidean(), B.locmat')
    mvt = pdf.(dist, distmat)
    mvt ./= sum(mvt, 2)
    MovementModel(mvt)
end

function MovementModel(B::Bathymetry,
    distance::Distributions.UnivariateDistribution,
    depth::Distributions.UnivariateDistribution)
    distmat = pairwise(Euclidean(), B.locmat')
    distmvt = pdf.(distance, distmat)
    dpthmvt = pdf.(depth, distmat)
    mvt = distmvt .* dpthmvt
    mvt ./= sum(mvt, 2)
    MovementModel(mvt)
end

"""
    eqdist(M::MovementModel)

Find the equilibrium distribution under a given movement model.
"""
function eqdist(M::MovementModel)
    eq = eigvecs(mvt')[:, 1]
    eq ./= sum(eq_dist)
    eq
end
function eqdist(M::MovementModel, B0::Real)
    eq = eqdist(M)
    B0 .* eqdist
end

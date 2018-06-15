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
    distmat = pairwise(Euclidean(), B.locs')
    distmvt = pdf.(distance, distmat)
    dpthmvt = pdf.(depth, B.bathy)
    mvt = distmvt .* dpthmvt'
    mvt ./= sum(mvt, 2)
    MovementModel(mvt)
end

"""
    eqdist(M::MovementModel)

Find the equilibrium distribution under a given movement model.
"""
function eqdist(M::MovementModel)
    eq = real(eigvecs(M.M')[:, 1])
    eq ./= sum(eq)
    PopState(eq)
end

(M::MovementModel)(P::PopState) = PopState(M.M' * P.P)

"""
    eqdist(M::MovementModel, B0::Real)

Finds the equilibrium spatial population distribution given a movement
operator `M` and total biomass `B0`. This is **slow**, because it
finds *all* of the eigenvalues and eigenvectors. I recommend using
`approx_eqdist` instead.
"""
function eqdist(M::MovementModel, B0::Real)
    eq = eqdist(M)
    PopState(B0 .* eq.P)
end

"""
    approx_eqdist(M::MovementModel, B0::Real)

Uses power iteration to find the spatial population distribution given
some movement operator `M` and total biomass `B0`. **Much** faster than
`eqdist` given reasonably large operators.
"""
function approx_eqdist(M::MovementModel, B0::Real)
    A = M.M'
    ## Tolerance taken from IterativeSolvers.jl
    tol = eps(real(eltype(A))) * size(A, 2) ^ 3
    d = size(A, 1)
    v0 = ones(d) / norm(ones(d))
    v1 = copy(v0)
    w = copy(v0)
    δv = 10 * tol * ones(d)
    while norm(δv, Inf) > tol
        w .= A * v0
        v1 .= w ./ norm(w)
        δv .= v0 .- v1
        v0 .= v1
    end
    PopState(B0 .* v1 ./ norm(v1, 1))
end

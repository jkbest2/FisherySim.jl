"""
    MovementModel{T<:Real}
        M::Matrix

Model for movement preferences.
"""
struct MovementModel{T<:Real, Ti<:Integer}
    M::Matrix{T}
    Ωsize::Tuple{Ti, Ti}
end

function MovementModel(domain::AbstractFisheryDomain,
                       habitat::AbstractMatrix,
                       pref_fn::Function,
                       dist_fn::Function)
    mvt_n = prod(domain.n)
    mvt = Matrix{eltype(habitat)}(undef, mvt_n, mvt_n)
    for jdx in 1:mvt_n, idx in 1:mvt_n
        mvt[idx, jdx] = pref_fn(habitat[idx]) *
            dist_fn(domain.distances[idx, jdx])
    end
    mvt ./= sum(mvt; dims = 1)
    MovementModel(mvt, domain.n)
end

function MovementModel(B::Bathymetry,
                       distance::Distributions.UnivariateDistribution)
    mvt = pdf.(dist, B.Ω.distances)
    mvt ./= sum(mvt, 2)
    MovementModel(Matrix(mvt'), size(B.Ω))
end

function MovementModel(B::Bathymetry,
                       distance::Distributions.UnivariateDistribution,
                       depth::Distributions.UnivariateDistribution)
    distmvt = pdf.(distance, B.Ω.distances)
    dpthmvt = pdf.(depth, vec(B.bathymetry))
    mvt = distmvt .* dpthmvt'
    mvt ./= sum(mvt; dims = 2)
    MovementModel(Matrix(mvt'), size(B.Ω))
end

function (M::MovementModel)(P::PopState)
    p = M.M * vecstate(P)
    PopState(reshape(p, M.Ωsize...))
end

"""
    eqdist(M::MovementModel)

Find the equilibrium distribution under a given movement model.
"""
function eqdist(M::MovementModel)
    n = M.Ωsize
    eq = real(eigvecs(M.M)[:, end])
    eq ./= sum(eq)
    PopState(reshape(eq, n...))
end

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

## TODO: Should stopping criterion rely on Rayleigh quotient giving an
## eigenvalue of one, rather than estimated eigenvector convergence?
## Also consider using ARPACK.jl (new, with FORTRAN deps) or
## IterativeSolvers.jl
"""
    approx_eqdist(M::MovementModel, B0::Real)

Uses power iteration to find the spatial population distribution given
some movement operator `M` and total biomass `B0`. **Much** faster than
`eqdist` given reasonably large operators.
"""
function approx_eqdist(M::MovementModel{T}) where T<:Real
    A = M.M
    n = M.Ωsize
    ## Tolerance taken from IterativeSolvers.jl, may need to be reduced.
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
    vnorm = v1 ./ norm(v1, 1)
    # @. v1 = v1 / norm(v1, 1)
    PopState(reshape(vnorm, n...))
end
function approx_eqdist(M::MovementModel{T}, B0::T) where T<:Real
    eq = approx_eqdist(M)
    PopState(eq.P .* B0)
end

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

Find the equilibrium distribution under a given movement model. Uses Arpack.jl
to find the eigenvector associated with the largest eigenvalue (which will be
one due to the row-standardization of the movement operator).
"""
function eqdist(M::MovementModel)
    n = M.Ωsize
    tol = 1e-10
    eq = zeros((0,))
    while true
        λ, ϕ = eigs(M.M; nev = 1, which = :LM, ritzvec = true,
                    tol = tol, maxiter = 1_000, v0 = eq)
        eq = vec(real(ϕ))
        if any(eq .< 0)
            tol *= 1e-1
        else
            break
        end
    end

    eq ./= sum(eq)
    PopState(reshape(eq, n...))
end

"""
    eqdist(M::MovementModel, B0::Real)

Finds the equilibrium spatial population distribution given a movement
operator `M` and total biomass `B0`.
"""
function eqdist(M::MovementModel, B0::Real)
    eq = eqdist(M)
    PopState(B0 .* eq.P)
end

struct MovementRate{F}
    dist::F

    function MovementRate(dist::F) where F<:Function
        new{F}(dist)
    end
end

(mr::MovementRate)(domain::AbstractFisheryDomain) = mr.dist.(domain.distances)

function MovementModel(habitat::Habitat{T, N} where T, habpref::HabitatPreference{F, N}, rate::MovementRate, domain::AbstractFisheryDomain) where {N, F}
    hpref = habpref(habitat; normalize = true)
    dist = rate(domain)

    movop = vec(hpref) .* dist
    movop ./= sum(movop; dims = 1)

    MovementModel(movop, size(domain))
end

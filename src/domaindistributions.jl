"""
    DomainDistribution{D<:Distribution, O<:AbstractFisheryDomain}
        dist::D
        domain::O

A distribution over a fishery domain. Methods for `rand` (without size) are
defined for univariate, multivariate, and `MatrixNormal` distributions, and
return realizations reshaped to the size of the fishery domain.
"""
struct DomainDistribution{D<:Distribution, O<:AbstractFisheryDomain}
    dist::D
    domain::O
end

## Construct a log-normal domain distribution where the μ vector provided
## defaults to the mean, but can be changed via the `s` argument (following the
## conventions of the `location` function from Distributions.jl).
function DomainDistribution(t::Type{<:Distributions.AbstractMvLogNormal},
                            μ::AbstractVector, Σ::AbstractMatrix,
                            domain::AbstractFisheryDomain,
                            s::Symbol = :mean)
    loc = location(t, s, μ, Σ)
    DomainDistribution(MvLogNormal(loc, Σ), domain)
end

function DomainDistribution(t::Type{<:Distributions.AbstractMvLogNormal},
                            μ::AbstractVector, Σ::AbstractPDMat,
                            domain::AbstractFisheryDomain,
                            s::Symbol = :mean)
    loc = location(t, s, μ, Σ.mat)
    DomainDistribution(MvLogNormal(loc, Σ), domain)
end

function DomainDistribution(t::Type{D},
                            μ::AbstractMatrix,
                            U::AbstractPDMat, V::AbstractPDMat,
                            domain::AbstractFisheryDomain,
                            s::Symbol = :mean) where {D<:MatrixLogNormal}
    loc = location(D, s, μ, U, V)
    normal = MatrixNormal(loc, U, V)
    DomainDistribution(t(normal), domain)
end

function rand(DD::DomainDistribution{D, O}) where {D<:UnivariateDistribution, O}
    rand(DD.dist, size(DD.domain)...)
end

function rand(DD::DomainDistribution{D, O}) where {D<:MultivariateDistribution, O}
    x = rand(DD.dist)
    reshape(x, size(DD.domain))
end

function rand(DD::DomainDistribution{D, O}) where {D<:MatrixDistribution, O}
    xs = rand(DD.dist)
    [reshape(x, size(DD.domain)) for x in eachcol(xs)]
    # reshape(xs, size(DD.domain)..., size(xs, 2))
end



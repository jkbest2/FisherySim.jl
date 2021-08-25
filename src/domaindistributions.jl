abstract type AbstractDomainDistribution end

"""
    DomainDistribution{D<:Distribution, O<:AbstractFisheryDomain}
        dist::D
        domain::O

A distribution over a fishery domain. Methods for `rand` (without size) are
defined for univariate, multivariate, and `MatrixNormal` distributions, and
return realizations reshaped to the size of the fishery domain.
"""
struct DomainDistribution{D, O<:AbstractFisheryDomain} <: AbstractDomainDistribution
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

function DomainDistribution(ddist::DomainDistribution, domain::AbstractFisheryDomain)
    DomainDistribution(ddist.dist, domain)
end

## Pass an instance, not a type (exceptions above)
function DomainDistribution(t::Type, domain::AbstractFisheryDomain)
    error("Passed $t as a type, need an instance e.g. $t()")
end

domain(dd::DomainDistribution) = dd.domain

function rand(DD::DomainDistribution{D, O}) where {D<:UnivariateDistribution, O}
    rand(DD.dist, size(DD.domain)...)
end

function rand(DD::DomainDistribution{D, O}) where {D<:MultivariateDistribution, O}
    x = rand(DD.dist)
    reshape(x, size(DD.domain))
end

function rand(dd::DomainDistribution{D, O}) where {D<:NeutralLandscapeMaker, O}
    rand(dd.dist, size(dd.domain))
end

# Blended distributions - combine multiple, weighted distributions over the domain
struct BlendedDomainDistribution{T} <: AbstractDomainDistribution
    ddists::Vector{AbstractDomainDistribution}
    wts::Vector{T}

    function BlendedDomainDistribution(ddists::Vector{AbstractDomainDistribution}, wts::Vector{T} = ones(length(ddists))) where {D<:AbstractDomainDistribution, T<:Real}
        length(ddists) == length(wts) || error("Must have weights for each distribution")
        all(==(domain(ddists[1])), domain.(ddists)) || error("All components must have the same domain")

        new{T}(ddists, wts)
    end
end

function BlendedDomainDistribution(dists::Vector, domain::AbstractFisheryDomain, wts::Vector{<:Real} = ones(length(dists)))
    ddist_vec = AbstractDomainDistribution[]
    for d in dists
        if d isa AbstractDomainDistribution
            push!(ddist_vec, d)
        else
            push!(ddist_vec, DomainDistribution(d, domain))
        end
    end

    BlendedDomainDistribution(ddist_vec, wts)
end

domain(bdd::BlendedDomainDistribution) = domain(bdd.ddists[1])

function rand(bddist::BlendedDomainDistribution)
    rvec = rand.(bddist.ddists)
    blend(rvec, bddist.wts)
end

# Classified domain distributions - classify
struct ClassifiedDomainDistribution{D, T, M} <: AbstractDomainDistribution
    ddist::D
    wts::Vector{T}
    mask::M

    function ClassifiedDomainDistribution(ddist::D, wts::Vector{T}, mask::M = nothing) where {D<:AbstractDomainDistribution, T<:Real, M}
        new{D, T, M}(ddist, wts, mask)
    end
end

function ClassifiedDomainDistribution(ddist::AbstractDomainDistribution, nclass::Integer, mask = nothing)
    wts = ones(nclass)
    ClassifiedDomainDistribution(ddist, wts, mask)
end

function ClassifiedDomainDistribution(dist, domain::AbstractFisheryDomain, classes, mask = nothing)
    ClassifiedDomainDistribution(DomainDistribution(dist, domain), classes, mask)
end

domain(cdd::ClassifiedDomainDistribution) = cdd.ddist.domain

function rand(cddist::ClassifiedDomainDistribution)
    res = rand(cddist.ddist)
    classify!(res, cddist.wts, cddist.mask)
end

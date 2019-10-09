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



function rand(DD::DomainDistribution{D, O}) where {D<:UnivariateDistribution, O}
    rand(DD.dist, size(DD.domain)...)
end

function rand(DD::DomainDistribution{D, O}) where {D<:MultivariateDistribution, O}
    x = rand(DD.dist)
    reshape(x, size(DD.domain))
end

function rand(DD::DomainDistribution{D, O}) where {D<:MatrixDistribution, O}
    xs = rand(DD.dist)
    reshape(xs, size(DD.domain)..., size(xs, 2))
end


struct MatrixLogNormal{T<:Real,M<:AbstractMatrix,C<:AbstractPDMat} <:
    Distribution{Matrixvariate,Continuous}
    normal::MatrixNormal{T, M, C}

    function MatrixLogNormal(normal::MatrixNormal{T, M, C}) where {T, M, C}
        new{T, M, C}(normal)
    end
end

function MatrixLogNormal(M, U, V)
    MatrixLogNormal(MatrixNormal(M, U, V))
end

rand(MLN::MatrixLogNormal) = exp.(rand(MLN.normal))

function location(::Type{D}, s::Symbol,
                  M::AbstractMatrix,
                  U::AbstractMatrix, V::AbstractMatrix) where {D<:MatrixLogNormal}
    loc = similar(M)
    @inbounds for jdx in size(loc, 2)
        M[:, jdx] .= location(MvLogNormal, s,
                              M[:, jdx],
                              diagm(V[jdx, jdx] .* diag(U)))
    end
    loc
end

function location(::Type{D}, s::Symbol,
                  M::AbstractMatrix,
                  U::AbstractPDMat, V::AbstractPDMat) where {D<:MatrixLogNormal}
    location(D, s, M, U.mat, V.mat)
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


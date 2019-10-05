abstract type AbstractCovarianceKernel end

check_covpars(varpar, lenpar) = varpar > 0 && lenpar > 0 ||
    throw(DomainError("Covariance parameters must be positive"))

"""
    ExpCov{T}
        σ²::T
        ϕ::T

The exponential covariance function with variance σ² and correlation distance
paramter ϕ, so that covariance at distance ``d`` is given by

```math
\\Cov(d) = \\sigma^2 \\exp\\left(-\\frac{d}{\\phi})
```

Note that this is equivalent to a Matérn covariance with smoothness parameter ν
= 1/2
"""
struct ExpCov{T} <: AbstractCovarianceKernel
    σ²::T
    ϕ::T

    function ExpCov(σ²::T, ϕ::T) where T
        check_covpars(σ², ϕ)
        new{T}(σ², ϕ)
    end
end

function ExpCov(σ², ϕ)
    σ², ϕ = promote(σ², ϕ)
    ExpCov(σ², ϕ)
end

function (K::ExpCov)(d)
    K.σ² * exp(-d / K.ϕ)
end

"""
    Matérn32Cov{T}
        σ²::T
        ϕ::T

The Matérn covariance function with fixed smoothness ν = 3/2. Covariance at distance ``d`` is given by

```math
\\Cov(d) = \\sigma^2 \\left(1 + \\frac{\\sqrt{3} d}{\\phi}\\right)
     \\exp\\left(-\\frac{\\sqrt{3} d}{\\phi}\\right)
```
"""
struct Matérn32Cov{T} <: AbstractCovarianceKernel
    σ²::T
    ϕ::T

    function Matérn32Cov(σ²::T, ϕ::T) where T
        check_covpars(σ², ϕ)
        new{T}(σ², ϕ)
    end
end

function Matérn32Cov(σ², ϕ)
    σ², ϕ = promote(σ², ϕ)
    Matérn32Cov(σ², ϕ)
end

# Add a non-Unicode name
const Matern32Cov = Matérn32Cov

function (K::Matérn32Cov)(d)
    K.σ² * (1 + √3 * d / K.ϕ) * exp(-√3 * d / K.ϕ)
end

function cov(K::AbstractCovarianceKernel, Ω::AbstractFisheryDomain)
    Σ = zeros(length(Ω), length(Ω))
    @inbounds for jdx in eachindex(Ω.locs), idx in eachindex(Ω.locs)
        idx > jdx && continue
        Σ[idx, jdx] = K(Ω.distances[idx, jdx])
    end
    PDMat(Symmetric(Σ, :U))
end


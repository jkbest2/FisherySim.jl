"""
    expcov(dist, σ, ϕ)

Exponential covariance kernel.
"""
function expcov(dist, σ, ϕ)
    σ * exp(-dist / ϕ)
end
"""
    sqexpcov(dist, σ, ϕ)

Squared exponential covariance kernel.
"""
function sqexpcov(dist, σ, ϕ)
    σ * exp(-dist^2 / 2 / ϕ)
end

"""
    BathymetryModel{T<:AbstractMvNormal}
        D::T
        locs::Array{Float64, 2}
        griddim::Vector{Float64}

A multivariate normal with defined R^2 locations for generating
bathymetry realizations
"""
struct BathymetryModel{Td<:AbstractMvNormal, Tω<:AbstractFisheryDomain} <: Any
    D::Td
    Ω::Tω
end
## Might be worth making a type system for spatial kernels...
function BathymetryModel(Ω::GriddedFisheryDomain{Tf, Ti},
                         μ::Vector{Tf},
                         kern::F) where {Tf<:Real, Ti, F<:Function}
    length(μ) == length(Ω) || throw(DimensionMismatch("Need mean for every grid cell."))
    D = MvNormal(μ, map_symm(kern, Ω.distances))
    BathymetryModel(D, Ω)
end
function BathymetryModel(Ω::GriddedFisheryDomain{Tf, Ti},
                         μ::Matrix{Tf},
                         kern::F) where {Tf<:Real, Ti, F<:Function}
    BathymetryModel(Ω, vec(μ), kern)
end

"""
    Bathymetry
        bathymetry::Matrix
        Ω::AbstractFisheryDomain

A realization of a `BathymetryModel`. Typically created using the
`rand(BM::BathyModel)` method.
"""
struct Bathymetry{T<:Real, Tω<:AbstractFisheryDomain} <: Any
    bathymetry::Matrix{T}
    Ω::Tω

    function Bathymetry(bathymetry::Matrix{T}, Ω::Tω) where {T<:Real, Tω<:AbstractFisheryDomain}
        all(size(Ω) .== size(bathymetry)) ||
            throw(DimensionMismatch("Bathymetry must conform to domain grid"))
        new{T, Tω}(bathymetry, Ω)
    end
end

function rand(rng::AbstractRNG, BM::BathymetryModel)
    bathy = rand(BM.D)
    bathy .-= minimum(bathy)
    Bathymetry(reshape(bathy, size(BM.Ω)...), BM.Ω)
end
rand(BM::BathymetryModel) = rand(Base.Random.GLOBAL_RNG, BM)

getindex(B::Bathymetry, i, j) = B.bathymetry[i, j]
getindex(B::Bathymetry, i) = B.bathymetry[i]
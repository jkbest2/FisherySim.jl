struct BathymetryModel{Td<:AbstractMvNormal, Tω<:AbstractFisheryDomain} <: Any
    D::Td
    Ω::Tω
end
## Might be worth making a type system for spatial kernels...
function BathymetryModel(Ω::GriddedFisheryDomain{Tf, Ti},
                         μ::AbstractVecOrMat{Tf},
                         kern::K) where {Tf<:Real, Ti, K<:AbstractCovarianceKernel}
    length(μ) == length(Ω) || throw(DimensionMismatch("Need mean for every grid cell."))
    D = MvNormal(vec(μ), cov(kern, Ω))
    BathymetryModel(D, Ω)
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

function rand(rng::Random.AbstractRNG, BM::BathymetryModel)
    bathy = rand(BM.D)
    bathy .-= minimum(bathy)
    Bathymetry(reshape(bathy, size(BM.Ω)...), BM.Ω)
end
rand(BM::BathymetryModel) = rand(Random.GLOBAL_RNG, BM)

getindex(B::Bathymetry, i, j) = B.bathymetry[i, j]
getindex(B::Bathymetry, i) = B.bathymetry[i]

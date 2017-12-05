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
struct BathymetryModel{T<:AbstractMvNormal} <: Any
    D::T
    locs::Array{Float64, 2}
    griddim::Vector{Float64}
end

"""
    Bathymetry
        bathy::Vector{Float64}
        locs::Array{Float64, 2}
        griddim::Vector{Float64}

A realization of a `BathymetryModel`. Typically created using the
`rand(BM::BathyModel)` method.
"""
struct Bathymetry <: Any
    bathy::Vector{Float64}
    locs::Array{Float64, 2}
    griddim::Vector{Float64}
end

function rand(BM::BathymetryModel)
    bathy = rand(BM.D)
    bathy .-= minimum(bathy)
    Bathymetry(bathy, BM.locs, BM.griddim)
end

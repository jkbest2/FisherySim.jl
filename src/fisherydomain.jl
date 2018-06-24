abstract type AbstractFisheryDomain <: Any end
abstract type DiscreteFisheryDomain <: AbstractFisheryDomain end 

struct GriddedFisheryDomain{T<:Real, Ti<:Integer} <: DiscreteFisheryDomain
    origin::Tuple{T, T}
    antipode::Tuple{T, T}
    n::Tuple{Ti, Ti}
    locs::Array{Tuple{T, T}, 2}
    distances::Array{T, 2}

    function GriddedFisheryDomain{T, Ti}(origin::Tuple{T, T},
                                  antipode::Tuple{T, T},
                                  n::Tuple{Ti, Ti},
                                  locs::Array{Tuple{T, T}, 2},
                                  distances::Array{T, 2}) where {T<:Real, Ti<:Integer}
        all(n .== size(locs)) ||
            throw(DimensionMismatch("Locations are wrong dimension."))
        all(prod(n) .== size(distances)) ||
            throw(DimensionMismatch("Distances are wrong dimension."))
        new(origin, antipode, n, locs, distances)
    end
end
"""
    GriddedFisheryDomain(origin::Tuple{T, T},
                         antipode::Tuple{T, T},
                         n::Tuple{Tn, Tn}) where {T<:Real, Tn<:Integer}

Construct a rectangular GriddedFisheryDomain with near and far corners 
origin and antipode respecitively, and n grid cells in each direction.
"""
function GriddedFisheryDomain(origin::Tuple{T, T},
                              antipode::Tuple{T, T},
                              n::Tuple{Tn, Tn}) where {T <: Real, Tn <: Integer}
    indices = Base.product(1:n[1], 1:n[2])
    ## Set grid up so that locations are at center of each grid cell
    domain_size = abs.(antipode .- origin)
    stepsize = domain_size ./ n
    offset = stepsize ./ 2
    locranges = range.(origin .+ offset, stepsize, n)
    locs = collect(Base.product(locranges...))

    distances = calculate_distances(locs)

    GriddedFisheryDomain{T, Tn}(origin, antipode, n, locs, distances)
end
function GriddedFisheryDomain(origin::Tuple{T, T},
                              antipode::Tuple{T, T},
                              n::Tn) where {T <: Real, Tn <: Integer}
    GriddedFisheryDomain(origin, antipode, (n, n))
end
size(Ω::GriddedFisheryDomain) = Ω.n
length(Ω::GriddedFisheryDomain) = prod(Ω.n)

function sample(rng, Ω::GriddedFisheryDomain, E::Integer; replace = true)
    N = prod(Ω.n)
    sample(rng, 1:N, E; replace = replace)
end
function sample(Ω::GriddedFisheryDomain, E::Integer; replace = true)
    sample(Base.Random.GLOBAL_RNG, Ω, E; replace = replace)
end
function sample(rng::AbstractRNG,
                Ω::GriddedFisheryDomain,
                w::StatsBase.AbstractWeights,
                E::Integer;
                replace = true)
    N = length(Ω)
    sample(rng, 1:N, w, E; replace = replace)
end
function sample(Ω::GriddedFisheryDomain, w::StatsBase.AbstractWeights, E::Integer; replace = true)
    sample(Base.Random.GLOBAL_RNG, Ω, w, E; replace = replace)
end

"""
    calculate_distances(locs::Array{Tuple{T, T}, 2})

Calcluate distances between points where locs is a matrix of tuples. Calculates
upper triangle of distances, then make symmetric, returning a Matrix{T, 2}.
"""
function calculate_distances(locs::Array{Tuple{T, T}, 2}) where T <: Real
    N = prod(size(locs))
    distances = Array{T, 2}(N, N)
    for j in 1:N
        for i in 1:(j - 1)
            distances[i, j] = hypot((locs[i] .- locs[j])...)
        end
        distances[j, j] = zero(T)
    end
    Matrix(Symmetric(distances))
end

"""
    map_symm(f::F, A::Array)

Map a function `f` over the elements of a symmetric matrix `A` by calculating
the upper-triangular elements and mirroring, returning the full matrix.
"""
function map_symm(f::F, A::Array{T, 2}) where {F<:Function, T}
    I, J = size(A)
    res = zeros(A)
    for j in 1:J, i in 1:j
        res[i, j] = f(A[i, j])
    end
    Matrix(Symmetric(res))
end

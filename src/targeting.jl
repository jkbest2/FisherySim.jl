abstract type AbstractTargetingBehavior <: Any end

"Vessels target fishing behavior at random within the fishery domain."
struct RandomTargeting <: AbstractTargetingBehavior end

function target(rng::Random.AbstractRNG,
                Ω::DiscreteFisheryDomain,
                t::RandomTargeting,
                E::Integer = 1)
    sample(rng, Ω, E; replace = true)
end

function target(Ω::DiscreteFisheryDomain,
                t::RandomTargeting,
                E::Integer = 1)
    target(Random.GLOBAL_RNG, Ω, t, E)
end

"Vessels target fixed locations in given order."
struct FixedTargeting{T} <: AbstractTargetingBehavior
    locations::Vector{T}
    ordered::Bool

    FixedTargeting(locs::Vector{T}) where T<:Integer = new{T}(locs, true)
end

function target(domain::AbstractFisheryDomain,
                t::FixedTargeting,
                E::Integer)
    E == length(t.locations) ||
        throw("Effort must equal number of fixed target locations.")
    t.locations
end

function target(domain::AbstractFisheryDomain,
                t::FixedTargeting)
    target(domain, t, length(t.locations))
end

## RNG is ignored, but allowed for compatibility with other targeting behaviors.
function target(rng::Random.AbstractRNG,
                domain::AbstractFisheryDomain,
                t::FixedTargeting,
                E::Integer)
    target(domain, t, length(t.locations))
end

"""
Vessels target according to some function or vector of preferences.

Accepts either a matrix of preference weights (don't necessarily need to add
to one) the same dimensions as the DiscreteFisheryDomain, or a function that
accepts a single argument (typically a Tuple or Array with two elements) of a
location and returns a preference weight.
"""
struct PreferentialTargeting{T} <: AbstractTargetingBehavior
    preference::T
end
function PreferentialTargeting(pref::Tp, Ω::DiscreteFisheryDomain) where {Tp <: AbstractArray}
    all(size(pref) .== size(Ω)) ||
        throw(DimensionMismatch("Preference array and domain dimensions must match."))
    PreferentialTargeting(pref)
end
function PreferentialTargeting(f::Function, Ω::DiscreteFisheryDomain)
    ## Quick @btime showed comprehension ~1/2 the time as `map`
    pref = [f(loc) for loc in Ω.locs]
    PreferentialTargeting(pref)
end

function target(rng::Random.AbstractRNG,
                Ω::DiscreteFisheryDomain,
                t::PreferentialTargeting{T},
                E::Integer = 1) where T <: AbstractArray
    sample(rng, Ω, Weights(vec(t.preference)), E; replace = true)
end
function target(Ω::DiscreteFisheryDomain,
                t::PreferentialTargeting{T},
                E::Integer = 1) where T <: AbstractArray
    target(Random.GLOBAL_RNG, Ω, t, E)
end

abstract type AbstractTargetingBehavior <: Any end

"Vessels target fishing behavior at random within the fishery domain."
struct RandomTargeting <: AbstractTargetingBehavior end

function target(rng::Random.AbstractRNG,
                t::RandomTargeting,
                Ω::DiscreteFisheryDomain)
    sample(rng, Ω)
end

function target(t::RandomTargeting,
                Ω::DiscreteFisheryDomain)
    target(Random.GLOBAL_RNG, t, Ω)
end

## FIXME: Should reimplement these (or at least FixedTargeting) as iterators
"""
Vessels target fixed locations in given order. This struct is mutable and has
state to track which survey station is next. State can be reset using the
`reset!` function.
"""
mutable struct FixedTargeting{T} <: AbstractTargetingBehavior
    locations::Vector{T}
    state::T

    FixedTargeting(locs::Vector{T}) where T<:Integer = new{T}(locs, one(T))
end

function target(t::FixedTargeting,
                domain::AbstractFisheryDomain)
    loc = t.locations[t.state]
    t.state += 1
    loc
end

## RNG is ignored, but allowed for compatibility with other targeting behaviors.
function target(rng::Random.AbstractRNG,
                t::FixedTargeting,
                domain::AbstractFisheryDomain)
    target(t, domain)
end

length(t::FixedTargeting) = length(t.locations)

"""
    reset!(FT::FixedTargeting{T}) where T
    reset!(Targeting::AbstractTargetingBehavior)

For `FixedTargeting` objects, resets `state` to one in place. For other
targeting types, does (and returns) nothing.
"""
reset!(FT::FixedTargeting{T}) where T = FT.state = one(T)
reset!(Targeting::AbstractTargetingBehavior) = nothing
reset!(target::AbstractTargetingBehavior, pop::PopState) = reset!(target)

abstract type AbstractPreferentialTargeting <: AbstractTargetingBehavior end

"""
Vessels target according to some function or vector of preferences.

Accepts either a matrix of preference weights (don't necessarily need to add
to one) the same dimensions as the DiscreteFisheryDomain, or a function that
accepts a single argument (typically a Tuple or Array with two elements) of a
location and returns a preference weight.
"""
struct PreferentialTargeting{T} <: AbstractPreferentialTargeting
    preference::T

    PreferentialTargeting(preference::W) where W<:Weights = new{W}(preference)
end

## Outer constructor fallback to construct Weights object for the inner constructor
PreferentialTargeting(preference) = PreferentialTargeting(weights(preference))

function PreferentialTargeting(pref::Tp, Ω::DiscreteFisheryDomain) where {Tp <: AbstractArray}
    all(size(pref) .== size(Ω)) ||
        throw(DimensionMismatch("Preference array and domain dimensions must match."))
    PreferentialTargeting(weights(pref))
end
function PreferentialTargeting(f::Function, Ω::DiscreteFisheryDomain)
    pref = f.(Ω.locs)
    PreferentialTargeting(weights(pref))
end

function target(rng::Random.AbstractRNG,
                t::AbstractPreferentialTargeting,
                Ω::DiscreteFisheryDomain)
    sample(rng, Ω, t.preference)
end
function target(t::AbstractPreferentialTargeting,
                Ω::DiscreteFisheryDomain)
    target(Random.GLOBAL_RNG, t, Ω)
end


mutable struct DynamicPreferentialTargeting{W, F} <: AbstractPreferentialTargeting
    preference::W
    pref_fn::F

    function DynamicPreferentialTargeting(preference::W, pref_fn::F) where
        {W<:Weights, F<:Function}
        new{W, F}(preference, pref_fn)
    end
end

function DynamicPreferentialTargeting(preference::AbstractArray,
                                      pref_fn::Function)
    DynamicPreferentialTargeting(weights(pref_fn(preference)), pref_fn)
end

function reset!(dynpref::DynamicPreferentialTargeting, pop::PopState)
    dynpref.preference = weights(dynpref.pref_fn(pop.P))
end


# Stratified random targeting; choose random sites within defined blocks
mutable struct StratifiedRandomTargeting{T} <: AbstractTargetingBehavior
    strata::Vector{Vector{T}}
    nsamples::Vector{T}
    state::T

    function StratifiedRandomTargeting(strata::Vector{Vector{T}}, nsamples::Vector{T}, state::T = one(T)) where T<:Integer
        any(nsamples .> 1) && error("Targeting doesn't work with >1 sample per strate yet.")
        new{T}(strata, nsamples, state)
    end
end

function StratifiedRandomTargeting(strata::Vector{Vector{T}}, n::T = 1) where T
    StratifiedRandomTargeting(strata, n * ones(T, length(strata)), one(T))
end

function StratifiedRandomTargeting(nblocks::Tuple, dom::GriddedFisheryDomain, n = 1)
    domsize = size(dom)
    dominds = LinearIndices(size(dom))
    blksize = domsize .÷ nblocks
    xinds, yinds = [1:b for b in blksize]

    strata = Vector{Vector{Int}}()
    for blk2 in 1:nblocks[2]
        colinds = xinds .+ blksize[2] * (blk2 - 1)
        for blk1 in 1:nblocks[1]
            rowinds = yinds .+ blksize[1] * (blk1 - 1)
            stratinds = vec(dominds[rowinds, colinds])
            push!(strata, stratinds)
        end
    end

    StratifiedRandomTargeting(strata, n)
end

reset!(stratrand::StratifiedRandomTargeting{T}) where T = stratrand.state = one(T)

function target(rng::Random.AbstractRNG,
                stratrand::StratifiedRandomTargeting,
                domain::AbstractFisheryDomain)
    loc = sample(rng, stratrand.strata[stratrand.state])
    stratrand.state += 1

    loc
end

function target(stratrand::StratifiedRandomTargeting,
                domain::AbstractFisheryDomain)
    target(Random.GLOBAL_RNG, stratrand, domain)
end

length(stratrand::StratifiedRandomTargeting) = sum(stratrand.nsamples)

abstract type AbstractTargetingBehavior <: Any end

"Vessels target fishing behavior at random within the fishery domain."
struct RandomTargeting <: AbstractTargetingBehavior end

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

function target(rng::AbstractRNG,
                Ω::DiscreteFisheryDomain,
                t::RandomTargeting,
                E::Integer = 1)
    sample(rng, Ω, E; replace = true)
end
function target(Ω::DiscreteFisheryDomain,
                t::RandomTargeting,
                E::Integer = 1)
    target(Base.Random.GLOBAL_RNG, Ω, t, E)
end
function target(rng::AbstractRNG,
                Ω::DiscreteFisheryDomain,
                t::PreferentialTargeting{T},
                E::Integer = 1) where T <: AbstractArray
    sample(rng, Ω, Weights(vec(t.preference)), E; replace = true)
end
function target(Ω::DiscreteFisheryDomain,
                t::PreferentialTargeting{T},
                E::Integer = 1) where T <: AbstractArray
    target(Base.Random.GLOBAL_RNG, Ω, t, E)
end

"""
Catchability describes a catchability surface for a vessel. Accepts either an
Array with the same dimensions as the DiscreteFisheryDomain or a function that
accepts a single argument (e.g. 2-element tuple or vector) of coordinates and
returns the catchability coefficient for that location. Note that at present
the function form is only used to form the array version under a particular
DiscreteFisheryDomain.
"""
struct Catchability{Tq}
    catchability::Tq
end
## If q is constant over the domain, just pass the value. Don't need to pass the
## DiscreteFisheryDomain, but it's allowed to keep the interface consistent.
Catchability(q::Tq) where Tq<:Real = Catchability{Tq}(q)
Catchability(q::Real, Ω::DiscreteFisheryDomain) = Catchability(q)
function Catchability(q::A, Ω::DiscreteFisheryDomain) where A <: AbstractArray
    all(size(q) .== size(Ω)) ||
        throw(DimensionMismatch("Catchability array and domain dimensions must match."))
    Catchability(q)
end
function Catchability(f::F, Ω::DiscreteFisheryDomain) where F <: Function
    q = [f(loc) for loc in Ω.locs]
    Catchability(q, Ω)
end


function getindex(Q::Catchability{Tq}, i, j) where Tq<:Real
    Q.catchability
end
function getindex(Q::Catchability{Tq}, i) where Tq<:Real
    Q.catchability
end
function getindex(Q::Catchability{Tq}, i, j) where Tq<:AbstractArray
    Q.catchability[i, j]
end
function getindex(Q::Catchability{Tq}, i) where Tq<:AbstractArray
    Q.catchability[i]
end

struct Vessel{Tt, Tq, Tf}
    target::Tt
    catchability::Tq
    ξ::Tf
    ϕ::Tf
    function Vessel(target::Tt,
                    catchability::Tq,
                    ξ::Tf,
                    ϕ::Tf) where {Tt<:AbstractTargetingBehavior,
                                  Tq<:Catchability,
                                  Tf<:Real}
        new{Tt, Tq, Tf}(target, catchability, ξ, ϕ)
    end
end

function Tweedie(μ::Tf, v::Vessel{Ta, Tq, Tf}) where {Ta, Tq, Tf<:Real}
     Tweedie(μ, v.ξ, v.ϕ)
end

## Vector of vessels needs to be abstract type for now; need to figure out small
## unions to type more concretely
struct Fleet{Tv, Te<:Integer}
    vessels::Vector{Tv}
    total_effort::Vector{Te}

    function Fleet(vessels::Vector{Tv}, total_effort::Vector{Te}) where {Tv<:Vessel, Te<:Integer}
        length(vessels) == length(total_effort) ||
            throw(DimensionMismatch("Must have an effort for each vessel"))
        new{Tv, Te}(vessels, total_effort)
    end
end

getindex(F::Fleet, i) = F.vessels[i]
vessels(F::Fleet) = F.vessels
length(F::Fleet) = length(F.vessels)

"""
    Catch
        time::Ti
        loc_idx::Ti
        coordinates::Tuple{Tf, Tf}
        effort::Tf
        c::Tf

Holds catch and effort data from a given fleet's fishing effort. Does
not account for within-season depletion.
"""
struct Catch{Tf, Ti} <: Any
    time::Ti
    loc_idx::Ti
    coordinates::Tuple{Tf, Tf}
    effort::Tf
    catch_biomass::Tf

    function Catch(time::Ti, loc_idx::Ti, coordinates::Tuple{Tf, Tf},
                   effort::Tf, catch_biomass::Tf) where {Tf<:Real, Ti<:Integer}
        new{Tf, Ti}(time, loc_idx, coordinates, effort, catch_biomass)
    end
end

"""
    fish!(P::PopState, V::Vessel, Ω::AbstractFisheryDomain, t::Integer)

Fish down a population state P. Returns a `Catch` object, mutates the
PopState in place. Limited to biomass available in a cell.
"""
function fish!(P::PopState,
               V::Vessel,
               Ω::AbstractFisheryDomain,
               t::Integer = 0)
    target_location = target(Ω, V.target)[1]
    μ = P.P[target_location] * V.catchability[target_location]
    if μ == 0
        catch_biomass = μ
    else
        catch_biomass = rand(Tweedie(μ, V.ξ, V.ϕ))
    end
    if catch_biomass > P.P[target_location]
        catch_biomass = P[target_location]
    end
    setindex!(P, P[target_location] - catch_biomass, target_location)
    C = Catch(t, target_location, Ω.locs[target_location],
              convert(typeof(catch_biomass), 1), catch_biomass)
end

"""
    fish!(P::PopState, F::Fleet, Ω::AbstractFisheryDomain, t::Integer)

Fish down a population P using vessels in fleet F and their associated efforts.
Vessels harvest in random order and cells are sequentially depleted.
"""
function fish!(P::PopState{Tf},
               F::Fleet,
               Ω::AbstractFisheryDomain,
               t::Ti = 0) where {Tf<:Real, Ti<:Integer}
    effort_vec = reduce(vcat, [repeat([vessel_idx], inner = [tot_eff]) for
                               (vessel_idx, tot_eff) in enumerate(F.total_effort)])
    catch_record = Vector{Catch{Tf, Ti}}()
    shuffle!(effort_vec)
    for idx in effort_vec
        c = fish!(P, F[idx], Ω, t)
        push!(catch_record, c)
    end
    catch_record
end

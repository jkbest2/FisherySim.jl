"""
    Vessel
        target::AbstractTargetingBehavior
        catchability::Catchability
        ξ::Real
        ϕ::Real

Describes the behavior of a fishing vessel, including its targeting strategy
(random or preferential), its catchability coefficient (which may vary across
space), and the parameters of the Tweedie distribution that govern its catch
distribution.
"""
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

## Vector of vessels needs to be abstract type for now; need to figure out small
## unions to type more concretely
struct Fleet{Tv, Te<:Integer}
    vessels::Vector{Tv}
    total_effort::Vector{Te}
    priority::Vector{Te}

    function Fleet(vessels::Vector{Tv},
                   total_effort::Vector{Te},
                   priority::Vector{Te}) where {Tv<:Vessel, Te<:Integer}
        length(vessels) == length(total_effort) == length(priority) ||
            throw(DimensionMismatch("Must have an effort and priority for each vessel"))
        priority = denserank(priority)
        new{Tv, Te}(vessels, total_effort, priority)
    end
end

function Fleet(vessels::Vector{Tv},
               total_effort::Vector{Te}) where {Tv<:Vessel, Te<:Integer}
    Fleet(vessels, total_effort, ones(eltype(total_effort), length(total_effort)))
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
    vessel_idx::Ti
    loc_idx::Ti
    coord1::Tf
    coord2::Tf
    effort::Tf
    catch_biomass::Tf

    function Catch(time::Ti, vessel_idx::Ti,
                   loc_idx::Ti, coord1::Tf, coord2::Tf,
                   effort::Tf, catch_biomass::Tf) where {Tf<:Real, Ti<:Integer}
        new{Tf, Ti}(time, vessel_idx, loc_idx, coord1, coord2, effort, catch_biomass)
    end
end
function Catch(time::Ti, vessel_idx::Ti,
               loc_idx::Ti, coordinates::Tuple{Tf, Tf},
               effort::Tf, catch_biomass::Tf) where {Tf<:Real, Ti<:Integer}
    Catch(time, vessel_idx, loc_idx, coordinates[1], coordinates[2], effort, catch_biomass)
end

"""
    fish!(P::PopState, V::Vessel, Ω::AbstractFisheryDomain, t::Integer)

Fish down a population state P. Returns a `Catch` object, mutates the
PopState in place. Limited to biomass available in a cell.
"""
function fish!(P::PopState,
               V::Vessel,
               Ω::AbstractFisheryDomain,
               t::Integer = 1,
               vessel_idx = 0)
    target_location = (l = target(V.target, Ω), t = t)
    μ = P.P[target_location.l] * V.catchability[target_location]
    if μ == 0
        catch_biomass = μ
    else
        catch_biomass = rand(CompoundPoissonGamma(μ, V.ξ, V.ϕ))
    end
    if catch_biomass > P.P[target_location.l]
        catch_biomass = P[target_location.l]
    end
    setindex!(P, P[target_location.l] - catch_biomass, target_location.l)
    C = Catch(t, vessel_idx, target_location.l, Ω.locs[target_location.l],
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
               t::Ti = 1) where {Tf<:Real, Ti<:Integer}
    effort_vec = order_effort(F)
    catch_record = Vector{Catch{Tf, Ti}}()
    shuffle!(effort_vec)
    for idx in effort_vec
        c = fish!(P, F[idx], Ω, t, idx)
        push!(catch_record, c)
    end
    # FIXME need a better way to deal with resetting stateful targeting each
    # year and dealing with changing targeting over time
    for vessel in vessels(F)
        reset!(vessel.target, P)
    end
    catch_record
end

"Return a vector with vessel indices ordered for fishing."
function order_effort(fleet::Fleet{Tv, Te}) where {Tv, Te}
    tot_eff = fleet.total_effort
    pri = fleet.priority
    vessel_idx = 1:length(fleet)
    eff_vvec = [fill(i, tot_eff[i]) for i in vessel_idx]

    eff = [Te[] for _ in unique(pri)]
    for (i, p) in enumerate(pri)
        append!(eff[p], eff_vvec[i])
    end
    shuffle!.(eff)
    reduce(vcat, eff)
end

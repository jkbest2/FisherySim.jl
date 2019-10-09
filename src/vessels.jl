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
    vessel_idx::Ti
    loc_idx::Ti
    coordinates::Tuple{Tf, Tf}
    effort::Tf
    catch_biomass::Tf

    function Catch(time::Ti, vessel_idx::Ti,
                   loc_idx::Ti, coordinates::Tuple{Tf, Tf},
                   effort::Tf, catch_biomass::Tf) where {Tf<:Real, Ti<:Integer}
        new{Tf, Ti}(time, vessel_idx, loc_idx, coordinates, effort, catch_biomass)
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
               t::Integer = 0,
               vessel_idx = 0)
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
    C = Catch(t, vessel_idx, target_location, Ω.locs[target_location],
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
        c = fish!(P, F[idx], Ω, t, idx)
        push!(catch_record, c)
    end
    catch_record
end

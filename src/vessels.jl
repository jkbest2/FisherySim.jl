abstract type Vessel <: Any end

"""
    SurveyVessel

Vessel in a survey fleet; samples at random within the region.
"""
struct SurveyVessel <: Vessel
    q::Float64
    Etot::Int
end
"""
    target(::Vessel, ::Popstate)

Choose set of fishing locations for a given vessel.
"""
function target(SV::SurveyVessel, P::PopState)
    sample(1:length(P.P), SV.Etot, replace = true)
    #rand(Multinomial(SV.Etot, length(P.P)))
end

"""
    FisheryVessel

Vessel that targets based on underlying population density.
"""
struct FisheryVessel <: Vessel
    q::Float64
    Etot::Int
end
function target(FV::FisheryVessel, P::PopState)
    sample(1:length(P.P), Weights(P.P), FV.Etot, replace = true)
    # rand(Multinomial(FV.Etot, P.P ./ sum(P.P)))
end

"""
    Catch

Holds catch and effort data from a given fleet's fishing effort. Does
not account for within-season depletion.
"""
struct Catch <: Any
    locs::Vector{Int}
    E::Vector{Float64}
    F::Vector{Float64}
end
function Catch(V::Vessel, P::PopState, σ::Float64)
    effort = target(V, P)
    ctch = zeros(Float64, effort)
    for (loc, eff) in enumerate(effort)
        ctch[loc] = rand(LogNormal(log(P.P[loc] * V.q * effort), σ))
    end
    Catch(effort, ctch)
end

# FIXME: add a Fleet type to deal with this better?
"""
    fish(P::PopState, Fleet::Vector{Vessel}, σ::Float64)

Fish the current population with the fleet in `VV`. Catches
are removed from the current population in random order over
the course of a season. Currently assumes constant effort at
each location that is fished (i.e. effort is 1 exactly).
"""
function fish(P::PopState, Fleet::Vector{Vessel}, σ::Float64)
    nv = length(Fleet)
    nloc = length(P.P)

    # Figure out which cells will be fished this season for each
    # vessel
    effort = target.(Fleet, P)

    ctch = Vector{Vector{Float64}}(nv)
    eff = Vector{Vector{Float64}}(nv)
    # Need to be able to index into each vessel for param values.
    v_idx = deepcopy(effort)
    for v in 1:nv
        fill!(v_idx[v], v)
        ctch[v] = zeros(Float64, nloc)
        eff[v] = zeros(Float64, nloc)
    end

    eff_vec, vvec = vcat(effort...), vcat(v_idx...)
    fish_order = randperm(length(vvec))

    Pnext = deepcopy(P)
    for trip in fish_order
        loc = eff_vec[trip]
        v = vvec[trip]
        if Pnext.P[loc] > 0
            c = rand(LogNormal(log(Pnext.P[loc] * Fleet[v].q), σ))
        else
            c = 0.0
        end
        if Pnext.P[loc] > c
            Pnext.P[loc] -= c
            ctch[v][loc] += c
        else
            ctch[v][loc] = Pnext.P[loc]
            Pnext.P[loc] = 0.0
        end
        eff[v][loc] += 1
    end

    Crec = Vector{Catch}(2)
    for v in 1:nv
        fished = eff[v] .> 0
        locvec = collect(1:length(eff[v]))[fished]
        Crec[v] = Catch(locvec, eff[v][fished], ctch[v][fished])
    end
    Pnext, Crec
end

struct CPUE
    locs::Vector{Int}
    cpue::Vector{Float64}
end
CPUE(C::Catch) = CPUE(C.locs, C.F ./ C.E)

# function +(C1::Catch, C2::Catch)
#     Catch(C1.E .+ C2.E, C1.F .+ C2.F)
# end

function -(P::PopState, C::Catch)
    PopState(P.P .- C.F)
end

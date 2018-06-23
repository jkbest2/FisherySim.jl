"""
    simulate(P::PopState,
             F::Fleet,
             M::MovementModel,
             PopDy::PopulationDynamicsModel,
             Ω::AbstractFisheryDomain,
             T::Int)

Simulate the population from P0 forward `T` years, by:

1. Remove individuals through fishing mortality using the `Fleet`,
2. Allow individuals to move according to their preferences,
3. Step the regional Schaefer model forward, distributing new
   individuals proportionally.
"""
function simulate(P0::PopState{Tf},
                  F::Fleet,
                  M::MovementModel,
                  PopDy::PopulationDynamicsModel,
                  Ω::AbstractFisheryDomain,
                  T::Ti) where {Tf<:Real, Ti<:Integer}
    P = copy(P0)
    Pvec = Vector{typeof(P)}()
    Cvec = Vector{Catch{Tf, Ti}}()
    for t in 1:T
        push!(Pvec, copy(P))
        c = fish!(P, F, Ω, t)
        append!(Cvec, c)
        P = M(P)                # TODO: Make these two steps in-place?
        P = step(PopDy, P)
    end
    Pvec, Cvec
end

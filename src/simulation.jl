"""
    simulate(P0::PopState,
             M::MovementModel,
             S::Schaefer,
             Fleet::Vector{Vessels},
             T::Int, σ::Float64)

Simulate the population from P0 forward `T` years, by:

1. Remove individuals through fishing mortality using the `Fleet`,
2. Allow individuals to move according to their preferences,
3. Step the regional Schaefer model forward, distributing new
   individuals proportionally.
"""
function simulate(P0::PopState,
                  M::MovementModel,
                  S::PopulationDynamicsModel,
                  Fleet::Vector{Vessel},
                  T::Int, σ::Float64)
    Pvec = Vector{PopState}(T + 1)
    Pvec[1] = P0
    Crecord = Vector{Vector{Catch}}(T)

    for yr in 2:(T + 1)
        Pf, Crecord[yr - 1] = fish(Pvec[yr - 1], Fleet, σ)
        Pm = M(Pf)
        Pvec[yr] = step(S, Pm)
    end

    Pvec, Crecord
end

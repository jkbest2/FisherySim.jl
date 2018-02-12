abstract type PopulationDynamicsModel end

"""
    PopState
        P::Vector{Float64}

Hold the population state at a given time.
"""
struct PopState <: Any
    P::Vector{Float64}
end

sum(P::PopState) = sum(P.P)

"""
    Schaefer
        r::Float64
        K::Float64

Parameters of a Schaefer population dynamics model.
"""
struct Schaefer <: PopulationDynamicsModel
    r::Float64
    K::Float64
end

"""
    step(S::Schaefer, p::PopState)

Step a population dynamics model forward by one time increment.
Population change is based on a region-wide carrying capacity,
but each cell steps forward individually.
"""
function step(S::Schaefer, P::PopState)
    Ptot = sum(P)
    Pnew = P.P .+  S.r .* P.P * (1 - Ptot / S.K)
    PopState(Pnew)
end

"""
    SchaeferStoch
        r::Float64
        K::Float64
        D::Distribution

Schaefer model with multiplicative process variation as
described by D.
"""
struct SchaeferStoch <: PopulationDynamicsModel
    r::Float64
    K::Float64
    D::Distribution
end

function step(S::SchaeferStoch, P::PopState)
    Pnew = step(Schaefer(S.r, S.K), P)
    PopState(rand(S.D, length(Pnew.P)) .* Pnew.P)
end

"""
    SchaeferKStoch
        r::Float64
        Kdist::Distibution

Schaefer model with random deviations of K
"""
struct SchaeferKStoch <: PopulationDynamicsModel
    r::Float64
    Kdist::Distribution
end

function step(S::SchaeferKStoch, P::PopState)
    K = rand(S.Kdist)
    step(Schaefer(S.r, K), P)
end

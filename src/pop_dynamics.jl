"""
    PopState
        P::Vector{Float64}

Hold the population state at a given time.
"""
struct PopState <: Any
    P::Vector{Float64}
end

"""
    Schaefer
        r::Float64
        K::Float64

Parameters of a Schaefer population dynamics model.
"""
struct Schaefer
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

sum(P::PopState) = sum(P.P)

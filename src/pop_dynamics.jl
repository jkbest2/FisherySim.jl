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
    step

Step a population dynamics model forward by one time increment. Here
assumes region-level dynamics and carrying capacity. New individuals
are allocated to locations based on original population.
"""
function step(S::Schaefer, P::PopState)
    Ptot = sum(P.P)
    Pnew = Ptot + S.r * (1 - Ptot / S.K)
    PopState(Pnew * P.P ./ Ptot)
end

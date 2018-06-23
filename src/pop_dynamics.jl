abstract type PopulationDynamicsModel end

"""
    PopState{T<:Real}
        P::Matrix{T}

Hold the population state at a given time.
"""
struct PopState{T<:Real} <: Any
    P::Matrix{T}
end

vecstate(P::PopState) = vec(P.P)
sum(P::PopState) = sum(P.P)

getindex(P::PopState, i) = P.P[i]
getindex(P::PopState, i, j) = P.P[i, j]
setindex!(P::PopState, x, i) = setindex!(P.P, x, i)
setindex!(P::PopState, x, i, j) = setindex!(P.P, x, i, j)

"""
    Schaefer{T<:Real}
        r::T
        K::T

Parameters of a Schaefer population dynamics model.
"""
struct Schaefer{T<:Real} <: PopulationDynamicsModel
    r::T
    K::T
end

## FIXME: This is a bad pun on `step`, which is used to get the step size of a
## `Range` object. Should come up with a better name.
"""
    step(S::Schaefer, p::PopState)

Step a population dynamics model forward by one time increment.
Population change is based on a region-wide carrying capacity,
but each cell steps forward individually.
"""
function step(S::Schaefer, P::PopState)
    Ptot = sum(P)
    Pnew = copy(P.P)
    real_r = S.r * (1 - Ptot / S.K)
    @. Pnew = P.P + real_r * P.P
    PopState(Pnew)
end

"""
    SchaeferStoch
        r::Float64
        K::Float64
        D::Distribution

Schaefer model with multiplicative process variation as
described by D, which must have a non-negative support.
"""
struct SchaeferStoch{T<:Real, Td<:Distribution} <: PopulationDynamicsModel
    r::T
    K::T
    D::Td

    function SchaeferStoch(r::T, K::T, D::Td) where {T<:Real, Td<:Distribution}
        minimum(support(D)) ≥ 0 ||
            throw(DomainError("Noise distribution must have non-negative support"))
        new{T, Td}(r, K, D)
    end
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

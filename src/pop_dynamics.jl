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
copy(P::PopState) = PopState(copy(P.P))

length(P::PopState) = length(P.P)
size(P::PopState) = size(P.P)
size(P::PopState, n::Integer) = size(P.P, n)

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

"""
    (S::Schaefer)(p::PopState)

Step a population state forward by one time increment. Population change is
based on a region-wide carrying capacity, but each cell steps forward
individually.
 """
function (S::Schaefer)(P::PopState)
    Ptot = sum(P)
    Pnew = copy(P)
    real_r = S.r * (1 - Ptot / S.K)
    # Here and below, can't use @. macro because then Pnew.P is lowered to
    # getfield.(Pnew, :P), which doesn't work (or make sense).
    Pnew.P .= P.P .+ real_r .* P.P
    Pnew
end

"""
    PellaTomlinson{T<:Real}
        r::T
        K::T
        m::T

A Pella-Tomlinson population dynamics model with growth rate `r`, carrying
capacity `K`, and shape parameter `m`. Population will change over time
following

```math
P_{t+1} = P_t + \\frac{r}{m-1} P_t \\left(1 - \\left(\\frac{\\sum P_t}{K}\\right)^{m-1}\\right)
```

Note that the Schaefer model is a special case of this model where ``m = 2``,
and the Fox model is a limiting case where ``m \\to 1``.
"""
struct PellaTomlinson{T<:Real} <: PopulationDynamicsModel
    r::T
    K::T
    m::T
end

"""
    (PT::PellaTomlinson)(P::PopState)

Step a population state forward by one time increment. Population change is
based on a region-wide carrying capacity but each cell steps forward
independently.
"""
function (PT::PellaTomlinson)(P::PopState)
    Ptot = sum(P)
    Pnew = copy(P)
    real_r = PT.r / (PT.m - 1) * (1 - (Ptot / PT.K) ^ (PT.m - 1))
    @. Pnew.P = P.P + real_r * P.P
    Pnew
end

"""
    StochasticProduction{M<:PopulationDynamicsModel, D<:Distribution}
        dynmod::M
        dist::D

Use `dynmod` population dynamics, and then apply multiplicative noise in the
resulting population.
"""
struct StochasticProduction{M<:PopulationDynamicsModel,
                            D<:Union{Distribution,DomainDistribution}} <:PopulationDynamicsModel
    dynmod::M
    dist::D
end

function (sp::StochasticProduction{M,D})(P::PopState) where {M,D<:DomainDistribution}
    Pnew = sp.dynmod(P)
    Pnew.P .= Pnew.P .* rand(sp.dist)
    Pnew
end

function (sp::StochasticProduction{M,D})(P::PopState) where {M,D<:UnivariateDistribution}
    Pnew = sp.dynmod(P)
    Pnew.P .= Pnew.P .* rand(sp.dist, size(Pnew)...)
    Pnew
end

function (sp::StochasticProduction{M,D})(P::PopState) where {M,D<:MultivariateDistribution}
    Pnew = sp.dynmod(P)
    Pnew.P .= Pnew.P .* reshape(rand(sp.dist), size(Pnew)...)
    Pnew
end


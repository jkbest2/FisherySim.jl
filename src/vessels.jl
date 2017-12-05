
abstract type Vessel <: Any end

struct SurveyVessel <: Vessel
    q::Float64
    Etot::Int
end
function target(SV::SurveyVessel, P::PopState)
    rand(Multinomial(length(P.P), FV.Etot))
end

struct FisheryVessel <: Vessel
    q::Float64
    Etot::Int
end
function target(FV::FisheryVessel, P::PopState)
    rand(Multinomial(FV.Etot, P.P ./ sum(P.P)))
end

struct Catch <: Any
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

struct CPUE
    cpue::Vector{Float64}
end
CPUE(C::Catch) = C.F ./ C.E

function +(C1::Catch, C2::Catch)
    Catch(C1.E .+ C2.E, C1.F .+ C2.F)
end

function -(P::PopState, C::Catch)
    PopState(P.P .- C.F)
end

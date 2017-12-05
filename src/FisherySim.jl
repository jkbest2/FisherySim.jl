module FisherySim

using Distributions
using StatsBase
using Distances

import Base: rand, +, -

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand

include("movement.jl")
export MovementModel, eqdist

include("pop_dynamics.jl")
export PopState, Schaefer, step

include("vessels.jl")
export Vessel, SurveyVessel, FisheryVessel, target, Catch, CPUE, +, -

end # module

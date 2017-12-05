module FisherySim

using Distributions
using StatsBase
using Distances

import Base: rand, +, -, step

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand

include("pop_dynamics.jl")
export PopState, Schaefer, step

include("movement.jl")
export MovementModel, eqdist

include("vessels.jl")
export Vessel, SurveyVessel, FisheryVessel, target, Catch, CPUE, +, -

end # module

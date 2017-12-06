module FisherySim

using Distributions
using StatsBase
using Distances

import Base: rand, +, -, step, sum

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand

include("pop_dynamics.jl")
export PopState, Schaefer, step, sum

include("movement.jl")
export MovementModel, eqdist

include("vessels.jl")
export Vessel, SurveyVessel, FisheryVessel, target, Catch, CPUE, +, -, fish

include("simulation.jl")
export simulate

end # module

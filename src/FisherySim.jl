module FisherySim

using Distributions
using Rtweedie
using StatsBase
using Distances
using IterativeSolvers

import Base: rand, +, -, step, sum

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand

include("pop_dynamics.jl")
export PopulationDynamicsModel, PopState, Schaefer, step, sum, SchaeferStoch, SchaeferKStoch

include("movement.jl")
export MovementModel, eqdist, approx_eqdist

include("vessels.jl")
export AbstractTargetingBehavio,
       RandomTargeting,
       PreferentialTargeting,
       target,
       Catchability,
       Vessel,
       Catch,
       logistic,
       CPUE,
       +, -,
       fish

include("simulation.jl")
export simulate

end # module

module FisherySim

using Distributions
using Rtweedie
using StatsBase
using Distances
using IterativeSolvers

import Base: rand, +, -, step, sum, getindex, size
import StatsBase: sample

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand

include("fisherydomain.jl")
export AbstractFisheryDomain,
       DiscreteFisheryDomain,
       GriddedFisheryDomain,
       size,
       sample
       ## Don't export:
       ## calculate_distances,
       ## map_symm

include("pop_dynamics.jl")
export PopulationDynamicsModel, PopState, Schaefer, step, sum, SchaeferStoch, SchaeferKStoch

include("movement.jl")
export MovementModel, eqdist, approx_eqdist

include("vessels.jl")
export AbstractTargetingBehavior,
       RandomTargeting,
       PreferentialTargeting,
       target,
       Catchability,
       Vessel,
       Catch,
       logistic,
       CPUE,
       +, -,
       fish,
       getindex

include("simulation.jl")
export simulate

end # module

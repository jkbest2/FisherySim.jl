module FisherySim

using Distributions
using Rtweedie
using StatsBase
using StatsFuns         # For `logit` and `logistic`; worth it?
using LinearAlgebra
using Random

import Base: rand, +, -, step, sum, getindex, setindex!, size, length, copy
import Rtweedie: Tweedie
import StatsBase: sample

include("fisherydomain.jl")
export AbstractFisheryDomain,
       DiscreteFisheryDomain,
       GriddedFisheryDomain,
       size,
       length,
       sample
       ## Don't export:
       ## calculate_distances,
       ## map_symm

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand

include("pop_dynamics.jl")
export PopulationDynamicsModel,
       PopState,
       Schaefer,
       vecstate,
       step,
       sum,
       SchaeferStoch,
       SchaeferKStoch,
       setindex!,
       copy

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
       fish!,
       getindex,
       Fleet,
       vessels

include("simulation.jl")
export simulate

end # module
